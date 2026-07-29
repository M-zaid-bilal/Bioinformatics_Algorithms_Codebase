import sys
import itertools
from Hamming_Distance import hamming_distance


def d_pattern_sequence(pattern, sequence):
    """
    Computes min Hamming distance between Pattern and any k-mer window in sequence.
    """
    k = len(pattern)
    min_dist = float('inf')
    
    for i in range(len(sequence) - k + 1):
        window = sequence[i:i+k]
        dist = hamming_distance(pattern, window)
        if dist < min_dist:
            min_dist = dist
            
    return min_dist


def d_pattern_dna(pattern, dna):
    """
    Computes sum of d(Pattern, sequence) across all DNA strings.
    """
    return sum(d_pattern_sequence(pattern, seq) for seq in dna)


def median_string(dna, k):
    """
    Finds the k-mer pattern that minimizes d(Pattern, Dna).
    """
    distance = float('inf')
    median = ""

    # Generate candidate k-mers in lexicographical order (AA...AA to TT...TT)
    for p_tuple in itertools.product('ACGT', repeat=k):
        pattern = "".join(p_tuple)
        current_d = d_pattern_dna(pattern, dna)
        
        if distance > current_d:
            distance = current_d
            median = pattern

    return median


def main():
    # File names (or pass via sys.argv if preferred)
    input_filename = "dataset_30304_9.txt"
    output_filename = "output_dataset_30304_9.txt"

    # Read dataset from file
    with open(input_filename, "r") as f:
        lines = [line.strip() for line in f if line.strip()]

    k = int(lines[0])
    dna = lines[1].split()

    # Solve Median String Problem
    result = median_string(dna, k)

    # Write output to file
    with open(output_filename, "w") as f:
        f.write(result + "\n")

    print(f"Result written to '{output_filename}': {result}")


if __name__ == "__main__":
    main()