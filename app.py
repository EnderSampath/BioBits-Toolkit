from flask import Flask, render_template, request
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from codon_table import CODON_TABLE
from codon_utils import (clean_sequence, split_into_codons, translate_codons, count_codon_usage, calculate_gc_content, plot_codon_usage, count_amino_acid_frequency, plot_amino_acid_frequency)
from tools import reverse_complement  
from tools import orf_finder
from tools.restriction_utils import simulate_digest_biopython
from Bio.Restriction import AllEnzymes
from tools import gc_skew

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        sequence = request.form.get('sequence', '')
    else:
        sequence = ''
    return render_template('index.html', sequence=sequence)

@app.route('/analyze', methods=['POST'])
def analyze():
    sequence = request.form.get('sequence', '').upper()

    if 'file' in request.files:
        file = request.files['file']
        if file.filename != '':
            sequence = file.read().decode('utf-8')

    graph_type = request.form.get('graph_type', 'codon')

    cleaned_seq = clean_sequence(sequence)
    codons = split_into_codons(cleaned_seq)
    protein = translate_codons(codons)
    gc_content = calculate_gc_content(cleaned_seq)

    # Select graph filename
    filename = 'codon_usage.png' if graph_type == 'codon' else 'aa_freq.png'
    filepath = os.path.join('static', filename)

    if graph_type == 'codon':
        codon_count = count_codon_usage(codons)
        plot_codon_usage(codon_count, save_path=filepath)
    else:
        aa_count = count_amino_acid_frequency(codons)

        # DEBUG prints to check what's happening
        print("🚨 Codons passed to AA count:", codons)
        print("🚨 Amino acid count:", aa_count)

        if not aa_count:
            print("⚠️ Warning: No valid amino acids found. Check the input sequence.")

        plot_amino_acid_frequency(aa_count, save_path=filepath)
    return render_template('results.html',
                           protein=protein,
                           gc_content=gc_content,
                           plot_filename=filename,
                           graph_type=graph_type)

@app.route('/reverse-complement', methods=['GET', 'POST'])
def reverse_complement_route():
    result = None
    plot_filename = None  # default

    if request.method == 'POST':
        sequence = request.form.get('sequence', '')
        cleaned = reverse_complement.clean_sequence(sequence)
        complement = reverse_complement.get_complement(cleaned)
        rev_complement = reverse_complement.get_reverse_complement(cleaned)
        gc = reverse_complement.calculate_gc_content(cleaned)

        result = {
            'original': cleaned,
            'complement': complement,
            'reverse_complement': rev_complement,
            'gc_content': gc
        }

        # ✅ ✅ ✅ Base composition bar chart code:
        import matplotlib.pyplot as plt
        import os

        base_counts = {
            'A': cleaned.count('A'),
            'T': cleaned.count('T'),
            'G': cleaned.count('G'),
            'C': cleaned.count('C')
        }

        fig, ax = plt.subplots()
        ax.bar(base_counts.keys(), base_counts.values(), color=['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'])
        ax.set_title('Base Composition')
        ax.set_ylabel('Count')
        plt.tight_layout()

        plot_filename = 'reverse_complement_composition.png'
        plot_path = os.path.join('static', plot_filename)
        plt.savefig(plot_path)
        plt.close()

    return render_template(
        'reverse_complement.html',
        result=result,
        plot_filename=plot_filename
    )

@app.route('/orf-finder', methods=['GET', 'POST'])
def orf_finder_route():
    results = []
    plot_filename = None

    if request.method == 'POST':
        sequence = request.form.get('sequence', '')
        results = orf_finder.find_orfs(sequence)

        if results:
            plot_path = os.path.join('static', 'orf_plot.png')
            orf_finder.plot_orfs(results, save_path=plot_path)
            plot_filename = 'orf_plot.png'

    return render_template("orf_finder.html", results=results, plot_filename=plot_filename)


from tools.restriction_utils import simulate_digest_biopython

@app.route("/restriction-enzyme", methods=["GET", "POST"])
def restriction():
    all_enzyme_names = sorted([enz.__name__ for enz in AllEnzymes])

    if request.method == "POST":
        sequence = request.form["sequence"]
        enzymes = request.form.getlist("enzymes")
        results, message = simulate_digest_biopython(sequence, enzymes)
        return render_template("restriction.html", enzymes=all_enzyme_names, results=results, message=message)

    return render_template("restriction.html", enzymes=all_enzyme_names)

@app.route('/gc-skew', methods=['GET', 'POST'])
def gc_skew_route():
    result = None
    if request.method == 'POST':
        sequence = request.form.get('sequence', '')
        cleaned_seq = gc_skew.clean_sequence(sequence)
        positions, skew_values = gc_skew.calculate_gc_skew(cleaned_seq)
        filename = 'gc_skew_plot.png'
        filepath = os.path.join('static', filename)
        gc_skew.plot_gc_skew(positions, skew_values, filepath)
        result = {
            'cleaned': cleaned_seq,
            'plot_filename': filename
        }
    return render_template('gc_skew.html', result=result)

from flask import send_from_directory
import os

@app.route('/resources')
def resources_index():
    return render_template('resources.html')

@app.route('/resources/download/<path:filename>')
def resources_download(filename):
    base_dir = os.path.join(os.path.dirname(__file__), 'resources', 'USABO')
    return send_from_directory(base_dir, filename, as_attachment=True)


# ✅ LAST lines of the file
if __name__ == '__main__':
    app.run(debug=True)
