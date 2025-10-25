def clean_sequence(seq):
    return ''.join(filter(lambda x: x in 'ATGC', seq.upper()))

def calculate_gc_skew(sequence, window_size=10):
    skew_values = []
    positions = []

    for i in range(0, len(sequence) - window_size + 1, window_size):
        window = sequence[i:i+window_size]
        g = window.count('G')
        c = window.count('C')

        if g + c == 0:
            skew = 0
        else:
            skew = (g - c) / (g + c)

        skew_values.append(skew)
        positions.append(i + window_size // 2)

    return positions, skew_values

def plot_gc_skew(positions, skew_values, save_path):
    import matplotlib.pyplot as plt

    print("Plotting GC Skew:")
    print("Positions:", positions)
    print("Skew values:", skew_values)

    plt.figure(figsize=(10, 4))
    plt.plot(positions, skew_values, color='purple', marker='o')  # ← Adds points to help visualize better
    plt.axhline(0, color='gray', linestyle='--')
    plt.title("GC Skew Across DNA Sequence")
    plt.xlabel("Position (bp)")
    plt.ylabel("GC Skew")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
