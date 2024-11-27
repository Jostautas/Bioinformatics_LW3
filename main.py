import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

def detect_quality_encoding(file_path):
    encodings = {
        "Sanger Phred+33": (33, 73),
        "Solexa Solexa+64": (59, 104),
        "Illumina 1.3+ Phred+64": (64, 104),
        "Illumina 1.5+ Phred+64": (3, 41),
        "Illumina 1.8+ Phred+33": (33, 74)
    }

    encoding_counts = {encoding: 0 for encoding in encodings}

    with open(file_path, 'r') as fastq_file:
        _ = fastq_file.readline()

        for quality_scores in fastq_file:
            for encoding, (min_val, max_val) in encodings.items():
                encoding_counts[encoding] += all(min_val <= ord(char) <= max_val for char in quality_scores.strip())

        detected_encoding = max(encoding_counts, key=encoding_counts.get)
        detected_range = encodings[detected_encoding]

    return detected_encoding, detected_range

def create_cg_proportion_histogram(file_path):
    proportions = []

    with open(file_path, 'r') as fastq_file:
        in_sequence = False
        sequence = ''
        for line in fastq_file:
            line = line.strip()
            if line.startswith('@'):
                in_sequence = True
                sequence = ''
            elif in_sequence:
                sequence += line
                in_sequence = False 

                cg_count = sequence.count('C') + sequence.count('G')
                total_nucleotides = len(sequence)
                proportion = (cg_count / total_nucleotides) * 100
                proportions.append((sequence, proportion))

    # Sort the proportions by C/G proportion in descending order
    proportions.sort(key=lambda x: x[1], reverse=True)

    return proportions

def find_top_sequences(file_path):
    sequences = []  # List to store the top sequences for each peak

    with open(file_path, 'r') as fastq_file:
        in_sequence = False
        sequence = []
        peak_start = 0
        peak_end = 0
        sequence_count = 0

        for line in fastq_file:
            line = line.strip()
            if line.startswith('@'):
                in_sequence = True
                sequence_count += 1
                sequence = [line]  # Start with the header
                peak_start = sequence_count
            elif in_sequence:
                sequence.append(line)
                if len(sequence) == 4:  # Once you have collected all four lines of a sequence
                    sequence_count += 1
                    sequences.append(''.join(sequence))
                    peak_end = sequence_count

                    # Check if this sequence is part of a new peak
                    if peak_end - peak_start == 4:
                        in_sequence = False

    return sequences

def save_sequences_to_fasta(file_name, sequences):
    with open(file_name, 'w') as fasta_file:
        for i, sequence in enumerate(sequences):
            fasta_file.write(f'>Sequence_{i + 1}\n')
            fasta_file.write(sequence + '\n')

file_path = 'reads_for_analysis.fastq'

# Detect the quality encoding
detected_encoding, detected_range = detect_quality_encoding(file_path)

if detected_encoding:
    print(f"Detected Quality Encoding: {detected_encoding}")
    print(f"Range for Detected Encoding: ({detected_range[0]}, {detected_range[1]})")
else:
    print("Could not determine the quality encoding.")

# Create the C/G proportion histogram
proportions = create_cg_proportion_histogram(file_path)

plt.figure(figsize=(8, 4))
plt.hist([p[1] for p in proportions], bins=np.linspace(0, 100, 300), edgecolor='k', alpha=0.6)
plt.xlabel('Proportion of C/G Nucleotides %')
plt.ylabel('Reads Frequency')
plt.title('Histogram of C/G Proportions')
plt.show()

# Find the top sequences
top_sequences = find_top_sequences(file_path)

# Separate the sequences into three peaks (assuming 3 peaks)
peak1_sequences = top_sequences[:5]
peak2_sequences = top_sequences[5:10]
peak3_sequences = top_sequences[10:15]

# Save sequences to FASTA files
save_sequences_to_fasta('peak_1_sequences.fasta', peak1_sequences)
save_sequences_to_fasta('peak_2_sequences.fasta', peak2_sequences)
save_sequences_to_fasta('peak_3_sequences.fasta', peak3_sequences)
