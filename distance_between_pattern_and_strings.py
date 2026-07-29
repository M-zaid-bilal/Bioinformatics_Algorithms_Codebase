import sys

try:
    from Hamming_Distance import hamming_distance
except ImportError:
    # Fallback hamming_distance implementation if Hamming_Distance.py is not in path
    def hamming_distance(p, q):
        """
        Calculates the Hamming distance between two equal-length strings p and q.
        """
        return sum(1 for a, b in zip(p, q) if a != b)


def distance_between_pattern_and_strings(pattern, dna):
    """
    Computes d(Pattern, Dna) = sum of min Hamming distances from Pattern to each string in Dna.
    
    Parameters:
    - pattern (str): Candidate k-mer pattern (e.g., 'AAA').
    - dna (list of str): Collection of DNA strings.
    
    Returns:
    - total_distance (int): Accumulated minimum Hamming distance across all sequences.
    """
    k = len(pattern)
    total_distance = 0

    for text in dna:
        min_hamming = float('inf')
        
        # Slide window across current DNA string
        for i in range(len(text) - k + 1):
            pattern_prime = text[i:i + k]
            dist = hamming_distance(pattern, pattern_prime)
            if min_hamming > dist:
                min_hamming = dist
                
        total_distance += min_hamming

    return total_distance


def solve_from_file(input_filename="dataset_dna_string_final.txt", output_filename="output_distance_pattern_dna_final.txt"):
    """
    Reads Pattern and Dna array from input_filename, computes d(Pattern, Dna),
    and writes the result to output_filename.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            print("Error: Input file must contain Pattern on line 1 and DNA strings on line 2.")
            return None

        pattern = lines[0]
        # DNA strings can be space-separated on line 2 or across multiple lines
        dna_line = lines[1]
        dna = dna_line.split() if " " in dna_line else lines[1:]

        result = distance_between_pattern_and_strings(pattern, dna)

        with open(output_filename, "w") as f:
            f.write(str(result) + "\n")

        print(f"\nCalculated Distance d({pattern}, Dna): {result}")
        print(f"Written to '{output_filename}'")
        return result

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None


def solve_from_string(raw_data_str):
    """
    Parses a raw multi-line string dataset (e.g., copied directly from Cogniterra).
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    pattern = lines[0]
    dna = lines[1].split() if " " in lines[1] else lines[1:]

    result = distance_between_pattern_and_strings(pattern, dna)
    print(f"Distance Result: {result}")
    return result


if __name__ == "__main__":
    filename = input("Enter input filename (or press Enter for default 'dataset_distance_pattern_dna.txt'): ").strip()
    if not filename:
        filename = "dataset_dna_string_final.txt"

    solve_from_file(filename)