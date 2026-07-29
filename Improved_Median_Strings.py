import itertools
import sys

try:
    from distance_between_pattern_and_strings import distance_between_pattern_and_strings
except ImportError:
    # Fallback hamming_distance implementation
    def hamming_distance(p, q):
        return sum(1 for a, b in zip(p, q) if a != b)

    # Fallback distance_between_pattern_and_strings implementation
    def distance_between_pattern_and_strings(pattern, dna):
        k = len(pattern)
        total_distance = 0
        for text in dna:
            min_hamming = float('inf')
            for i in range(len(text) - k + 1):
                pattern_prime = text[i:i + k]
                dist = hamming_distance(pattern, pattern_prime)
                if min_hamming > dist:
                    min_hamming = dist
            total_distance += min_hamming
        return total_distance


def all_strings(k):
    """
    Generates all possible DNA strings of length k over the alphabet {A, C, G, T}.
    Returns a list of 4^k strings.
    """
    nucleotides = ['A', 'C', 'G', 'T']
    return [''.join(p) for p in itertools.product(nucleotides, repeat=k)]


def median_string(dna, k):
    """
    Finds a k-mer Pattern that minimizes d(Pattern, Dna) across all possible k-mers.
    
    Parameters:
    - dna (list of str): Collection of DNA sequences.
    - k (int): Length of the k-mer pattern to search for.
    
    Returns:
    - median (str): A k-mer pattern with the minimal total distance to Dna.
    """
    distance = float('inf')
    median = ""
    patterns = all_strings(k)

    for pattern in patterns:
        current_dist = distance_between_pattern_and_strings(pattern, dna)
        # Strict inequality ensures we return the first k-mer in case of ties
        if distance > current_dist:
            distance = current_dist
            median = pattern

    return median


def solve_from_file(input_filename="dataset_median_string.txt", output_filename="output_median_string.txt"):
    """
    Reads k and DNA sequences from input_filename, runs median_string(),
    and writes the resulting median k-mer to output_filename.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            print("Error: Input file must contain k on line 1 and DNA sequences on line 2.")
            return None

        k = int(lines[0])
        dna_line = lines[1]
        dna = dna_line.split() if " " in dna_line else lines[1:]

        result = median_string(dna, k)

        with open(output_filename, "w") as f:
            f.write(result + "\n")

        print(f"\nMedian String (k={k}): {result}")
        print(f"Written to '{output_filename}'")
        return result

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None


def solve_from_string(raw_data_str):
    """
    Parses a raw multi-line string dataset (e.g. copied from Cogniterra).
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    k = int(lines[0])
    dna = lines[1].split() if " " in lines[1] else lines[1:]

    result = median_string(dna, k)
    print(f"Median String Result: {result}")
    return result


if __name__ == "__main__":
    filename = input("Enter input filename (or press Enter for default 'dataset_median_string.txt'): ").strip()
    if not filename:
        filename = "dataset_median_string.txt"

    solve_from_file(filename)