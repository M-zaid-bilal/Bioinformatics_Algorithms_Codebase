import random
import sys
import time

# Global map for fast nucleotide lookup to optimize hot loops
NUC_MAP = {'A': 0, 'C': 1, 'G': 2, 'T': 3}


def create_profile_with_pseudocounts(motifs):
    """
    Creates a 4 x k Profile matrix from Motifs using Laplace's Rule of Succession (pseudocounts of 1).
    Divided by (len(motifs) + 4) to normalize counts to probabilities.
    """
    k = len(motifs[0])
    num_motifs = len(motifs)
    
    # Initialize count matrix with 1 for pseudocounts
    counts = [[1 for _ in range(k)] for _ in range(4)]
    
    # Accumulate nucleotide frequencies across motifs
    for motif in motifs:
        for j, nuc in enumerate(motif):
            counts[NUC_MAP[nuc]][j] += 1
            
    # Normalize counts to probabilities
    total_col_count = num_motifs + 4
    profile = [[counts[r][c] / total_col_count for c in range(k)] for r in range(4)]
    return profile


def pr_kmer_profile(kmer, profile):
    """
    Computes Pr(kmer | Profile) by multiplying position-specific nucleotide probabilities.
    """
    prob = 1.0
    for i, nuc in enumerate(kmer):
        prob *= profile[NUC_MAP[nuc]][i]
    return prob


def profile_randomly_generated_kmer(text, k, profile):
    """
    Generates a k-mer from text biased by probability distribution Pr(kmer | Profile).
    """
    num_kmers = len(text) - k + 1
    kmers = [text[i:i + k] for i in range(num_kmers)]
    
    # Calculate probability for each k-mer window
    probs = [pr_kmer_profile(kmer, profile) for kmer in kmers]
    total_prob = sum(probs)
    
    # Handle edge case where all probabilities are 0 (fallback to uniform distribution)
    if total_prob == 0:
        return random.choice(kmers)
    
    # Normalize probabilities
    norm_probs = [p / total_prob for p in probs]
    
    # Biased random selection based on weights
    selected_kmer = random.choices(kmers, weights=norm_probs, k=1)[0]
    return selected_kmer


def score_motifs(motifs):
    """
    Calculates total score of motifs (sum of non-consensus nucleotide counts per column).
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


def gibbs_sampler(dna, k, t, N):
    """
    Executes a single run of GibbsSampler for N outer iterations.
    
    Parameters:
    - dna: List of t DNA sequences.
    - k: Length of motif.
    - t: Number of DNA sequences.
    - N: Number of iterations per run.
    """
    # 1. Randomly select initial k-mers in each sequence
    motifs = []
    for seq in dna:
        start_idx = random.randint(0, len(seq) - k)
        motifs.append(seq[start_idx:start_idx + k])
        
    best_motifs = list(motifs)

    # 2. Iterate N times
    for _ in range(N):
        # Pick a random sequence index i to remove
        i = random.randint(0, t - 1)
        
        # Build Profile from all motifs except Motifi
        motifs_except_i = [motifs[j] for j in range(t) if j != i]
        profile = create_profile_with_pseudocounts(motifs_except_i)
        
        # Sample new Motifi from sequence Dna[i] using Profile distribution
        motifs[i] = profile_randomly_generated_kmer(dna[i], k, profile)
        
        # Update best_motifs if strict improvement
        if score_motifs(motifs) < score_motifs(best_motifs):
            best_motifs = list(motifs)

    return best_motifs


def run_gibbs_sampler(dna, k, t, N, num_starts=20):
    """
    Runs GibbsSampler multiple times (default: 20 starts) with random initializations
    to avoid local minima traps.
    """
    start_time = time.time()
    best_motifs = gibbs_sampler(dna, k, t, N)
    best_score = score_motifs(best_motifs)

    print(f"\n[+] Executing GibbsSampler ({num_starts} runs, N={N} iterations each)...")

    for r in range(1, num_starts):
        current_motifs = gibbs_sampler(dna, k, t, N)
        current_score = score_motifs(current_motifs)
        
        if current_score < best_score:
            best_score = current_score
            best_motifs = current_motifs

    elapsed = time.time() - start_time
    print(f" -> Completed {num_starts} runs in {elapsed:.2f}s | Best Score: {best_score}")

    return best_motifs


def solve_from_file(input_filename="dataset_gibbs_sampler.txt", output_filename="output_gibbs_sampler.txt", N=1000, num_starts=20):
    """
    Reads k, t, N and DNA sequences from input_filename, runs GibbsSampler,
    and writes space-separated results to output_filename.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        first_line = lines[0].split()
        k = int(first_line[0])
        t = int(first_line[1])
        if len(first_line) > 2:
            N = int(first_line[2])

        dna = []
        for line in lines[1:]:
            dna.extend(line.split())

        results = run_gibbs_sampler(dna, k, t, N, num_starts=num_starts)

        with open(output_filename, "w") as f:
            f.write(" ".join(results) + "\n")

        print(f"\nGibbsSampler Results:")
        print(" ".join(results))
        print(f"Score: {score_motifs(results)}")
        print(f"\nWritten to '{output_filename}'")
        return results

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None


def solve_from_string(raw_data_str, num_starts=20):
    """
    Parses a raw multi-line dataset string (e.g., copied from Cogniterra/Coursera).
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    first_line = lines[0].split()
    k = int(first_line[0])
    t = int(first_line[1])
    N = int(first_line[2]) if len(first_line) > 2 else 1000

    dna = []
    for line in lines[1:]:
        dna.extend(line.split())

    results = run_gibbs_sampler(dna, k, t, N, num_starts=num_starts)
    print(f"\nGibbsSampler Results: {' '.join(results)}")
    print(f"Score: {score_motifs(results)}")
    return results


if __name__ == "__main__":
    filename = input("Enter input filename (or press Enter for default 'dataset_gibbs_sampler.txt'): ").strip()
    if not filename:
        filename = "dataset_gibbs_sampler.txt"

    solve_from_file(filename, num_starts=20)