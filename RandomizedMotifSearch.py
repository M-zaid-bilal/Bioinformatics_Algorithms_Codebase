import random
import sys
import time

# Global map for fast nucleotide lookup to eliminate dictionary re-creation overhead in hot loops
NUC_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}

# Import profile_most_probable_kmer from profile_most_probable_kmer.py or profile_most_probable_k_mer.py
try:
    from profile_most_probable_kmer import profile_most_probable_kmer
except ImportError:
    try:
        from profile_most_probable_k_mer import profile_most_probable_kmer
    except ImportError:
        def pr_kmer_profile(kmer, profile):
            prob = 1.0
            for i, nuc in enumerate(kmer):
                prob *= profile[NUC_MAP[nuc]][i]
            return prob

        def profile_most_probable_kmer(text, k, profile):
            max_prob = -1.0
            best_kmer = text[:k]
            for i in range(len(text) - k + 1):
                window = text[i:i+k]
                prob = pr_kmer_profile(window, profile)
                if prob > max_prob:
                    max_prob = prob
                    best_kmer = window
            return best_kmer


def create_profile_with_pseudocounts(motifs):
    """
    Creates a 4 x k Profile matrix from Motifs using Laplace's Rule of Succession (pseudocounts of 1).
    """
    k = len(motifs[0])
    num_motifs = len(motifs)
    
    # Count matrix initialized with 1 for pseudocounts
    counts = [[1 for _ in range(k)] for _ in range(4)]
    
    for motif in motifs:
        for j, nuc in enumerate(motif):
            counts[NUC_MAP[nuc]][j] += 1
            
    # Normalize counts to probabilities
    total_col_count = num_motifs + 4
    profile = [[counts[r][c] / total_col_count for c in range(k)] for r in range(4)]
    return profile


def score_motifs(motifs):
    """
    Calculates total score of motifs (sum of non-consensus counts across columns).
    """
    k = len(motifs[0])
    num_motifs = len(motifs)
    score = 0

    for j in range(k):
        col_counts = [0, 0, 0, 0]
        for motif in motifs:
            col_counts[NUC_MAP[motif[j]]] += 1
        max_freq = max(col_counts)
        score += (num_motifs - max_freq)

    return score


def motifs_from_profile(profile, dna, k):
    """
    Generates a collection of Profile-most probable k-mers from each sequence in DNA.
    """
    return [profile_most_probable_kmer(seq, k, profile) for seq in dna]


def randomized_motif_search(dna, k, t):
    """
    Executes a single run of RandomizedMotifSearch until score stops improving.
    """
    # 1. Randomly select k-mers Motifs in each string from DNA
    motifs = []
    for seq in dna:
        start_idx = random.randint(0, len(seq) - k)
        motifs.append(seq[start_idx:start_idx + k])
    
    best_motifs = list(motifs)

    # 2. Iterate continuously until score stops improving
    while True:
        profile = create_profile_with_pseudocounts(motifs)
        motifs = motifs_from_profile(profile, dna, k)
        
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = list(motifs)
        else:
            return best_motifs


def run_randomized_motif_search(dna, k, t, num_iterations=1000):
    """
    Runs RandomizedMotifSearch N times (default: 1000) to avoid local minima,
    logging progress every 100 iterations to prevent perception of hangs.
    """
    start_time = time.time()
    best_motifs = randomized_motif_search(dna, k, t)
    best_score = score_motifs(best_motifs)

    print(f"\n[+] Starting {num_iterations} iterations of RandomizedMotifSearch...")

    for i in range(1, num_iterations):
        current_motifs = randomized_motif_search(dna, k, t)
        current_score = score_motifs(current_motifs)
        
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs

        # Print progress status every 100 iterations
        if (i + 1) % 100 == 0 or (i + 1) == num_iterations:
            elapsed = time.time() - start_time
            print(f" -> Progress: {i + 1}/{num_iterations} iterations | Current Best Score: {best_score} | Elapsed: {elapsed:.2f}s")

    return best_motifs


def solve_from_file(input_filename="dataset_randomized_motif.txt", output_filename="output_randomized_motif.txt", iterations=1000):
    """
    Reads k, t, and DNA sequences from an input file, runs RandomizedMotifSearch,
    and writes space-separated motifs to an output file.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        first_line = lines[0].split()
        k = int(first_line[0])
        t = int(first_line[1])

        dna = []
        for line in lines[1:]:
            dna.extend(line.split())

        results = run_randomized_motif_search(dna, k, t, num_iterations=iterations)

        with open(output_filename, "w") as f:
            f.write(" ".join(results) + "\n")

        print(f"\nRandomizedMotifSearch Results ({iterations} iterations):")
        print(" ".join(results))
        print(f"Score: {score_motifs(results)}")
        print(f"\nWritten to '{output_filename}'")
        return results

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None


def solve_from_string(raw_data_str, iterations=1000):
    """
    Parses a raw multi-line string dataset and runs RandomizedMotifSearch.
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    first_line = lines[0].split()
    k = int(first_line[0])
    t = int(first_line[1])

    dna = []
    for line in lines[1:]:
        dna.extend(line.split())

    results = run_randomized_motif_search(dna, k, t, num_iterations=iterations)
    print(f"\nRandomizedMotifSearch Result: {' '.join(results)}")
    print(f"Score: {score_motifs(results)}")
    return results


if __name__ == "__main__":
    filename = input("Enter input filename (or press Enter for default 'dataset_randomized_motif.txt'): ").strip()
    if not filename:
        filename = "dataset_rms.txt"

    solve_from_file(filename, iterations=1000)