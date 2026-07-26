# Importing functions from your existing modules
from Frequent_words_efficient import maxMap
from Freequent_words_with_mismatches import (
    approximate_pattern_count,
    approximate_pattern_count_from_file,
    frequent_words_with_mismatches,
    frequent_words_with_mismatches_from_file,
    neighbors,
)
from Reverse_compliment import reverse_complement


def frequent_words_with_mismatches_and_rc(text: str, k: int, d: int) -> str:
    """Finds all most frequent k-mers (with up to d mismatches and reverse complements)."""
    freq_map = {}
    n = len(text)

    # Populate frequency table considering both original and reverse complement neighborhoods
    for i in range(n - k + 1):
        pattern = text[i : i + k]
        pattern_rc = reverse_complement(pattern)

        # 1. Process neighborhood of original pattern
        pattern_neighbors = neighbors(pattern, d)
        for neighbor in pattern_neighbors:
            freq_map[neighbor] = freq_map.get(neighbor, 0) + 1

        # 2. Process neighborhood of reverse complement pattern
        rc_neighbors = neighbors(pattern_rc, d)
        for neighbor in rc_neighbors:
            freq_map[neighbor] = freq_map.get(neighbor, 0) + 1

    # Find maximum count and collect top k-mers
    max_count = maxMap(freq_map)
    frequent_patterns = [
        pattern for pattern, count in freq_map.items() if count == max_count
    ]

    return " ".join(frequent_patterns)


def frequent_words_with_mismatches_and_rc_from_file(
    input_filepath: str, output_filepath: str = None
) -> str:
    """Reads Text from line 1, and k, d from line 2 of the input file.

    Handles space-delimited k and d robustly.
    """
    with open(input_filepath, "r") as file:
        # Strip lines and filter out empty lines to avoid indexing issues
        lines = [
            line.strip() for line in file.read().splitlines() if line.strip()
        ]

    text = lines[0]

    # Split on any whitespace to handle space-separated integers "4 1"
    k_str, d_str = lines[1].split()
    k = int(k_str)
    d = int(d_str)

    result = frequent_words_with_mismatches_and_rc(text, k, d)

    if output_filepath:
        with open(output_filepath, "w") as out_file:
            out_file.write(result)
        print(f"Results successfully written to {output_filepath}")

    return result


# ==========================================
# Execution Example
# ==========================================
if __name__ == "__main__":
    # Test reading directly from file with the 2-line structure
    filepath = "dataset_30278_10.txt"
    try:
        file_result = frequent_words_with_mismatches_and_rc_from_file(filepath)
        print(f"Frequent Words with Mismatches & RC from '{filepath}':")
        print(file_result)
    except FileNotFoundError:
        print(f"File '{filepath}' not found.")