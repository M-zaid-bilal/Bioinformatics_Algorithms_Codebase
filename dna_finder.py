"""
dnaa_finder.py
Imports standalone local python modules to locate candidate DnaA boxes
from a raw genome text file and generates visualization plots.
"""
import json
import matplotlib.pyplot as plt
# Direct imports from individual standalone files in the same directory
from Frequent_Words_with_mismatches_and_reverse_complement import frequent_words_with_mismatches_and_rc
from Reverse_compliment import reverse_complement
from Hamming_Distance import hamming_distance
from Skew_Count import calculate_skew, minimum_skew


def read_sequence_from_txt(filepath: str) -> str:
    """Reads a raw sequence string from a .txt file, stripping whitespace/newlines."""
    with open(filepath, 'r') as file:
        return file.read().replace('\n', '').replace('\r', '').strip().upper()


def plot_skew_diagram(
    genome: str,
    ori_center: int,
    start_window: int,
    end_window: int,
    output_filename: str = "gc_skew_plot.png"
):
    """
    Calculates GC Skew across the entire genome, plots the line graph,
    highlights the minimum skew position (oriC), and saves the chart to a file.
    """
    # 1. Calculate GC Skew array using imported utility
    skew_values = calculate_skew(genome)

    # 2. Setup figure
    plt.figure(figsize=(12, 6))
    plt.plot(skew_values, label="GC Skew", color="#1f77b4", linewidth=1.5)

    # 3. Highlight oriC and window region
    plt.axvline(
        x=ori_center,
        color="red",
        linestyle="--",
        linewidth=2,
        label=f"oriC Center ({ori_center})"
    )
    plt.axvspan(
        start_window,
        end_window,
        color="orange",
        alpha=0.3,
        label=f"Search Window [{start_window}:{end_window}]"
    )

    # 4. Styling graph
    plt.title("Genome GC Skew Diagram & Estimated oriC Region", fontsize=14, fontweight="bold")
    plt.xlabel("Genome Position (bp)", fontsize=12)
    plt.ylabel("Skew Value (G - C)", fontsize=12)
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend(loc="upper right")
    plt.tight_layout()

    # 5. Save chart
    plt.savefig(output_filename, dpi=300)
    plt.close()
    print(f"[+] GC Skew plot saved to '{output_filename}'")


def find_dna_a_box(
    filepath: str,
    k: int = 9,
    d: int = 1,
    window_size: int = 500,
    plot_skew: bool = True
) -> dict:
    """
    Finds candidate DnaA boxes by running skew analysis and k-mer search.

    :param filepath: Path to the genome .txt file.
    :param k: Length of the DnaA box (default 9).
    :param d: Allowed Hamming distance mismatches (default 1).
    :param window_size: Search window size around the minimum skew origin.
    :param plot_skew: Boolean toggle to generate GC skew visualization graph.
    :return: Dictionary with candidate DnaA box motifs and their positions.
    """
    # 1. Load sequence
    genome = read_sequence_from_txt(filepath)

    # 2. Find oriC region via minimum skew
    min_skew_indices = minimum_skew(genome)
    ori_center = min_skew_indices[0]

    # 3. Crop window around minimum skew point
    start = max(0, ori_center - (window_size // 2))
    end = min(len(genome), ori_center + (window_size // 2))
    ori_window = genome[start:end]

    # 4. Find most frequent k-mers with mismatches + reverse complements
    candidate_motifs = frequent_words_with_mismatches_and_rc(ori_window, k, d).split()

    # 5. Map motif occurrences back to genome indices
    motif_locations = {}
    for motif in candidate_motifs:
        rc_motif = reverse_complement(motif)
        matches = []
        
        for i in range(len(ori_window) - k + 1):
            pattern = ori_window[i:i + k]
            if hamming_distance(pattern, motif) <= d or hamming_distance(pattern, rc_motif) <= d:
                matches.append(start + i)
                
        motif_locations[motif] = matches

    # 6. Plot GC Skew Diagram
    if plot_skew:
        plot_skew_diagram(genome, ori_center, start, end)

    return {
        "genome_length": len(genome),
        "ori_center": ori_center,
        "window_bounds": (start, end),
        "candidate_motifs": candidate_motifs,
        "motif_locations": motif_locations
    }


if __name__ == "__main__":
    target_file = "Salmonella_enterica.txt"

    print(f"[*] Attacking {target_file}...")
    
    try:
        results = find_dna_a_box(target_file, k=9, d=1, window_size=500, plot_skew=True)

        # Output to Terminal
        print("\n=== DnaA Box Analysis Results ===")
        print(f"Genome Length: {results['genome_length']} bp")
        print(f"Estimated oriC Center: {results['ori_center']}")
        print(f"Window Analyzed: {results['window_bounds'][0]} to {results['window_bounds'][1]}")
        print(f"Candidate DnaA Boxes: {results['candidate_motifs']}")
        
        print("\nLocations:")
        for motif, positions in results["motif_locations"].items():
            print(f"  > Motif '{motif}': {positions}")

        # Output to JSON File
        with open("dnaa_results.json", "w") as out:
            json.dump(results, out, indent=4)
        print("\n[+] Results saved to 'dnaa_results.json'")

    except FileNotFoundError:
        print(f"[!] Error: Target file '{target_file}' not found. Make sure it's in the same folder!")