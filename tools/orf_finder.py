import matplotlib.pyplot as plt

def clean_sequence(seq):
    return ''.join(filter(lambda x: x in 'ATGC', seq.upper()))

def find_orfs(sequence, min_length=30):
    sequence = clean_sequence(sequence)
    orfs = []
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    for frame in range(3):
        i = frame
        while i < len(sequence) - 2:
            codon = sequence[i:i+3]
            if codon == start_codon:
                for j in range(i + 3, len(sequence) - 2, 3):
                    stop = sequence[j:j+3]
                    if stop in stop_codons:
                        length = j + 3 - i
                        if length >= min_length:
                            orfs.append({
                                "start": i,
                                "end": j + 3,
                                "length": length,
                                "frame": frame,
                                "sequence": sequence[i:j+3]
                            })
                        break
                i = j + 3  # Skip past stop codon
            else:
                i += 3  # Move to next codon
    return orfs

def plot_orfs(orfs, save_path='static/orf_plot.png'):
    plt.figure(figsize=(10, 3))
    if not orfs:
        plt.text(0.5, 0.5, "No ORFs Found", fontsize=14, ha='center', va='center')
        plt.axis('off')
    else:
        y_ticks = []
        labels = []
        for idx, orf in enumerate(orfs):
            plt.barh(idx, orf['length'], left=orf['start'], height=0.6, color='skyblue', edgecolor='black')
            y_ticks.append(idx)
            labels.append(f"Frame {orf['frame']}")
        plt.yticks(y_ticks, labels)
        plt.xlabel("Nucleotide Position")
        plt.title("ORF Visualization")
        plt.grid(axis='x', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
