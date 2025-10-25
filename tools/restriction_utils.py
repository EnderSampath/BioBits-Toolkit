from Bio.Restriction import AllEnzymes, RestrictionBatch, Analysis
from Bio.Seq import Seq

def simulate_digest_biopython(sequence_str, enzyme_list):
    sequence = Seq(sequence_str.upper().replace(" ", "").replace("\n", ""))
    enzyme_class_map = {enz.__name__: enz for enz in AllEnzymes}

    try:
        selected_classes = [enzyme_class_map[name] for name in enzyme_list]
    except KeyError as e:
        return {}, f"Error: Unknown enzyme '{e.args[0]}'"

    batch = RestrictionBatch(selected_classes)
    cut_positions = batch.search(sequence)  # returns {enzyme_obj: [cut positions]}
    
    fragment_output = {}

    for enzyme_obj, cuts in cut_positions.items():
        enzyme_name = enzyme_obj.__name__

        # Sort and calculate fragments
        sorted_cuts = sorted(cuts)
        positions = [0] + sorted_cuts + [len(sequence)]
        fragment_lengths = [positions[i+1] - positions[i] for i in range(len(positions) - 1)]

        fragment_output[enzyme_name] = {
            "cut_sites": sorted_cuts,
            "fragment_lengths": fragment_lengths
        }

    return fragment_output, f"Successfully simulated with {len(selected_classes)} enzyme(s)."
