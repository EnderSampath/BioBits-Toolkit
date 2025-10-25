import matplotlib.pyplot as plt
from codon_table import CODON_TABLE
from collections import defaultdict

# Amino acid abbreviations and full names
AA_ABBREV = {
    'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
    'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
    'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
    'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val',
    'X': '???', '_': 'Stop'
}

AA_FULL = {
    'A': 'Alanine', 'R': 'Arginine', 'N': 'Asparagine', 'D': 'Aspartic acid', 'C': 'Cysteine',
    'Q': 'Glutamine', 'E': 'Glutamic acid', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'L': 'Leucine', 'K': 'Lysine', 'M': 'Methionine', 'F': 'Phenylalanine', 'P': 'Proline',
    'S': 'Serine', 'T': 'Threonine', 'W': 'Tryptophan', 'Y': 'Tyrosine', 'V': 'Valine',
    'X': 'Unknown', '_': 'Stop'
}

# === Utility Functions ===

def count_amino_acid_frequency(codons):
    from collections import defaultdict
    aa_count = defaultdict(int)
    for codon in codons:
        if len(codon) == 3 and codon in CODON_TABLE:
            aa = CODON_TABLE[codon]
            aa_count[aa] += 1
    return dict(aa_count)  # Convert defaultdict to normal dict to avoid Jinja issues

def plot_amino_acid_frequency(aa_count, save_path=None):
    labels = list(aa_count.keys())
    values = [aa_count[k] for k in labels]

    plt.figure(figsize=(12, 6))
    bars = plt.bar(labels, values, color='lightgreen', edgecolor='black')
    for i, v in enumerate(values):
        plt.text(i, v + 0.1, str(v), ha='center', fontsize=8)

    plt.title("Amino Acid Frequency")
    plt.xlabel("Amino Acids")
    plt.ylabel("Frequency")
    plt.tight_layout()

    if save_path:
        plt.savefig(save_path)
        print(f"✅ Amino acid graph saved to {save_path}")
    else:
        plt.show()  # Show the plot only if not saving

    plt.close()  # Close the plot to free memory

def load_sequence_from_file(filename):
    """Loads DNA sequence from a text file"""
    try:
        with open(filename, 'r') as file:
            return file.read()
    except FileNotFoundError:
        print("❌ File not found. Falling back to manual input.")
        return input("Enter a DNA sequence (ATGC only): ")

def clean_sequence(seq):
    """Cleans the DNA sequence: removes whitespace and converts to uppercase"""
    seq = seq.upper().replace(" ", "").replace("\n", "")
    return seq[:len(seq) - (len(seq) % 3)]  # trim to multiple of 3

def split_into_codons(seq):
    """Splits cleaned DNA sequence into codons (triplets)"""
    return [seq[i:i+3] for i in range(0, len(seq), 3)]

def translate_codons(codons):
    """Translates a list of codons into a protein sequence"""
    return ''.join([CODON_TABLE.get(codon, 'X') for codon in codons])

def count_codon_usage(codons):
    """Counts how frequently each codon appears"""
    codon_count = defaultdict(int)
    for codon in codons:
        codon_count[codon] += 1
    return codon_count

def print_codon_usage_table(codon_count):
    """Prints codon frequency table"""
    print("\n📊 Codon Usage:")
    for codon in sorted(codon_count):
        print(f"{codon}: {codon_count[codon]}")

def get_codon_label(codon, style='codon'):
    aa = CODON_TABLE.get(codon, 'X')
    if style == 'codon':
        return codon
    elif style == 'abbrev':
        return f"{codon} ({AA_ABBREV.get(aa, aa)})"
    elif style == 'full':
        return f"{codon} ({AA_FULL.get(aa, aa)})"
    else:
        return codon  # default fallback

def plot_codon_usage(codon_count, label_style='codon',save_path=None):
    """Plots codon usage bar graph with color highlights and chosen label style"""
    sorted_codons = sorted(codon_count.keys())
    frequencies = [codon_count[codon] for codon in sorted_codons]

    labels = [get_codon_label(codon, label_style) for codon in sorted_codons]

    colors = []
    for codon in sorted_codons:
        if codon == 'ATG':
            colors.append('green')  # Start codon
        elif codon in ['TAA', 'TAG', 'TGA']:
            colors.append('red')    # Stop codons
        else:
            colors.append('skyblue')

    plt.figure(figsize=(16, 7))
    plt.bar(labels, frequencies, color=colors, edgecolor='black')

    # Add frequency labels on top of bars
    for i, freq in enumerate(frequencies):
        plt.text(i, freq + 0.1, str(freq), ha='center', va='bottom', fontsize=8)

    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xlabel("Codons", fontsize=12)
    plt.ylabel("Frequency", fontsize=12)
    plt.title("Codon Usage Frequency", fontsize=14)
    plt.xticks(rotation=90, fontsize=10)
    plt.yticks(fontsize=10)
    plt.tight_layout()
    if save_path:
        plt.savefig(save_path)
        print(f"✅ Graph saved as {save_path}")
    plt.show()

def calculate_gc_content(seq):
    """Calculates GC content percentage of the DNA sequence"""
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    if total == 0:
        return 0
    return round(((g+c)/total)*100, 2)
# === Main Program ===

def main():
    print("🧬 Codon Usage & Protein Translator")

    # Ask how to input DNA sequence
    use_file = input("Do you want to load the DNA sequence from a file? (yes/no): ").strip().lower()
    if use_file in ['yes', 'y']:
        filename = input("Enter the filename (e.g. sequence.txt): ").strip()
        seq = load_sequence_from_file(filename)
    else:
        seq = input("Enter a DNA sequence (ATGC only): ")

    cleaned_seq = clean_sequence(seq)

    gc_content = calculate_gc_content(cleaned_seq)
    print(f"\n🧬 GC Content: {gc_content}%")
    
    codons = split_into_codons(cleaned_seq)
    protein = translate_codons(codons)
    codon_count = count_codon_usage(codons)

    print("\n🧪 Protein Sequence:")
    print(protein.upper())

    print_codon_usage_table(codon_count)

    # Ask user label style choice
    print("\nChoose label style for codon usage graph x-axis:")
    print("1 - Codon only (e.g., ATG)")
    print("2 - Codon + Amino Acid abbreviation (e.g., ATG (Met))")
    print("3 - Codon + Full Amino Acid name (e.g., ATG (Methionine))")
    label_choice = input("Enter 1, 2, or 3: ").strip()
    style_map = {'1': 'codon', '2': 'abbrev', '3': 'full'}
    label_style = style_map.get(label_choice, 'codon')

    save_graph = input("Do you want to save the graph as an image file? (yes/no): ").strip().lower()
    if save_graph in ['yes', 'y']:
        filename = input("Enter filename to save as (e.g., codon_usage.png): ").strip()
    else:
        filename = None

    plot_codon_usage(codon_count, label_style=label_style, save_path=filename)

    # Count amino acid frequency
    aa_count = count_amino_acid_frequency(codons)
   
    # Ask user if they want to save/show the amino acid frequency graph
    save_aa_graph = input("Do you want to save the amino acid frequency graph? (yes/no): ").strip().lower()
    if save_aa_graph in ['yes', 'y']:
        aa_filename = input("Enter filename to save amino acid graph (e.g. aa_frequency.png): ").strip()
    else:
        aa_filename = None

    plot_amino_acid_frequency(aa_count, save_path=aa_filename)



if __name__ == "__main__":
    main()

