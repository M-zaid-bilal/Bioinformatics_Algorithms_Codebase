import sys
from profile_most_probable_k_mer import profile_most_probable_kmer, pr_kmer_profile


def create_profile(motifs, k):
    """
    Creates a 4 x k profile matrix from a list of motif strings.
    Row 0: 'A', Row 1: 'C', Row 2: 'G', Row 3: 'T'.
    """
    t_count = len(motifs)
    profile = [[0.0] * k for _ in range(4)]
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for motif in motifs:
        for j, nuc in enumerate(motif):
            profile[nuc_map[nuc]][j] += 1.0

    # Convert counts to probabilities
    for r in range(4):
        for c in range(k):
            profile[r][c] /= t_count

    return profile


def score_motifs(motifs):
    """
    Calculates Score(Motifs) by summing the number of non-consensus 
    nucleotides in each column across all motifs.
    Lower score indicates a more conserved motif set.
    """
    k = len(motifs[0])
    t = len(motifs)
    total_score = 0
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

    for col in range(k):
        counts = [0] * 4
        for motif in motifs:
            counts[nuc_map[motif[col]]] += 1
        max_count = max(counts)
        total_score += (t - max_count)

    return total_score


def greedy_motif_search(dna, k, t):
    """
    Implements the standard GreedyMotifSearch algorithm (without pseudocounts).
    
    Input:
      - dna: List of t DNA sequence strings
      - k: Integer k-mer length
      - t: Number of sequences in dna
      
    Output:
      - BestMotifs: List of k-mers (one from each sequence in dna) minimizing Score(Motifs)
    """
    # BestMotifs initialized to first k-mer of each sequence
    best_motifs = [seq[:k] for seq in dna]
    best_score = score_motifs(best_motifs)

    # First sequence length
    n = len(dna[0])

    # Iterate over every k-mer in the first string of DNA
    for i in range(n - k + 1):
        motif_1 = dna[0][i:i+k]
        current_motifs = [motif_1]

        # Iteratively find Profile-most probable k-mer for remaining sequences
        for j in range(1, t):
            # Form profile from motifs formed so far
            profile = create_profile(current_motifs, k)
            
            # Find Profile-most probable k-mer in the j-th string
            motif_j = profile_most_probable_kmer(dna[j], k, profile)
            current_motifs.append(motif_j)

        current_score = score_motifs(current_motifs)

        # Strict inequality ensures tie-breaking selects earlier motifs
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs

    return best_motifs


# ... existing code ...
def solve_from_file(input_filename="dataset_gms.txt", output_filename="output_gms_2.txt"):
    """
    Reads k, t, and DNA sequences from an input file, runs greedy_motif_search,
    and writes space-separated motifs to an output file.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        first_line = lines[0].split()
        k = int(first_line[0])
        t = int(first_line[1])

        dna = lines[1:]
        # Handles case where DNA sequences are space-separated on second line
        if len(dna) == 1:
            dna = dna[0].split()

        results = greedy_motif_search(dna, k, t)

        with open(output_filename, "w") as f:
            f.write(" ".join(results) + "\n")

        print("GreedyMotifSearch Results:")
        print(" ".join(results))
        print(f"\nWritten to '{output_filename}'")
        return results

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None

if __name__ == "__main__":
    # Interactive file prompt
    filename = input("Enter input filename (e.g., dataset_gms.txt): ").strip()
    if not filename:
        filename = "dataset_gms.txt"

    solve_from_file(filename)