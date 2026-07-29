from profile_most_probable_k_mer import profile_most_probable_kmer
from greedy_motif_search import score_motifs


def create_profile_with_pseudocounts(motifs):
    """
    Creates a 4 x k Profile matrix from Motifs using Laplace's Rule of Succession.
    Adds pseudocount of 1 for A, C, G, T at every position.
    Divided by (len(motifs) + 4) to convert counts to probabilities.
    """
    k = len(motifs[0])
    num_motifs = len(motifs)
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    
    # Initialize count matrix with pseudocounts of 1
    counts = [[1 for _ in range(k)] for _ in range(4)]
    
    # Add observed counts from motifs
    for motif in motifs:
        for j, nuc in enumerate(motif):
            counts[nuc_map[nuc]][j] += 1
            
    # Normalize by total occurrences (num_motifs + 4 pseudocounts)
    total_col_count = num_motifs + 4
    profile = [[counts[r][c] / total_col_count for c in range(k)] for r in range(4)]
    
    return profile


def greedy_motif_search_with_pseudocounts(dna, k, t):
    """
    GreedyMotifSearch algorithm incorporating Laplace's Rule of Succession (pseudocounts).
    """
    # Form initial BestMotifs using the first k-mer in each DNA string
    best_motifs = [seq[:k] for seq in dna]

    # Loop through each candidate k-mer in the first DNA sequence
    for i in range(len(dna[0]) - k + 1):
        motif1 = dna[0][i:i+k]
        current_motifs = [motif1]

        for j in range(1, t):
            # Apply Laplace's Rule of Succession to form profile from current_motifs
            profile = create_profile_with_pseudocounts(current_motifs)
            
            # Find Profile-most probable k-mer in the j-th DNA string
            next_motif = profile_most_probable_kmer(dna[j], k, profile)
            current_motifs.append(next_motif)

        # Update best_motifs if strict improvement in score
        if score_motifs(current_motifs) < score_motifs(best_motifs):
            best_motifs = current_motifs

    return best_motifs


def solve_from_file(input_filename="dataset_greedy_ps.txt", output_filename="output_greedy_pseudocounts.txt"):
    """
    Reads k, t, and DNA sequences from an input file, runs greedy_motif_search_with_pseudocounts,
    and writes space-separated motifs to an output file.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        first_line = lines[0].split()
        k = int(first_line[0])
        t = int(first_line[1])

        dna = lines[1:]
        if len(dna) == 1:
            dna = dna[0].split()

        results = greedy_motif_search_with_pseudocounts(dna, k, t)

        with open(output_filename, "w") as f:
            f.write(" ".join(results) + "\n")

        print("\nGreedyMotifSearch (with Pseudocounts) Results:")
        print(" ".join(results))
        print(f"\nWritten to '{output_filename}'")
        return results

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None


if __name__ == "__main__":
    filename = input("Enter input filename (or press Enter for default 'dataset_greedy_ps.txt'): ").strip()
    if not filename:
        filename = "dataset_greedy_ps.txt"

    solve_from_file(filename)