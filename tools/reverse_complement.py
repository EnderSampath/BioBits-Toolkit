# tools/reverse_complement.py

def clean_sequence(seq):
    return seq.upper().replace(" ", "").replace("\n", "")

def get_complement(seq):
    complement_map = str.maketrans("ATGC", "TACG")
    return seq.translate(complement_map)

def get_reverse_complement(seq):
    complement = get_complement(seq)
    return complement[::-1]

def calculate_gc_content(seq):
    g = seq.count('G')
    c = seq.count('C')
    total = len(seq)
    if total == 0:
        return 0
    return round(((g + c) / total) * 100, 2)

import matplotlib.pyplot as plt
import os

def plot_base_composition(seq, save_path):
    counts = {
        'A': seq.count('A'),
        'T': seq.count('T'),
        'G': seq.count('G'),
        'C': seq.count('C')
    }

    labels = list(counts.keys())
    values = list(counts.values())
    colors = ['#a3c9a8', '#ffb7b2', '#89cff0', '#d1b3c4']

    plt.figure(figsize=(6, 4))
    plt.bar(labels, values, color=colors, edgecolor='black')
    plt.xlabel("Nucleotide")
    plt.ylabel("Count")
    plt.title("Base Composition of Reverse Complement")
    plt.tight_layout()
    plt.savefig(save_path)
    plt.close()
