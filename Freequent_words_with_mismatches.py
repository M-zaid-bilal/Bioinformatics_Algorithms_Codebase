# Importing functions from your existing modules
from Frequent_words_efficient import maxMap
from Hamming_Distance import hamming_distance


def approximate_pattern_count(text: str, pattern: str, d: int) -> int:
    """Counts occurrences of Pattern in Text with at most d mismatches."""
    count = 0
    pattern_len = len(pattern)

    for i in range(len(text) - pattern_len + 1):
        substring = text[i : i + pattern_len]
        if hamming_distance(pattern, substring) <= d:
            count += 1

    return count


def approximate_pattern_count_from_file(
    input_filepath: str, output_filepath: str = None
) -> int:
    """Reads Pattern (line 1), Text (line 2), and d (line 3) from input_filepath,

    calculates approximate pattern count, and optionally writes to output_filepath.
    """
    with open(input_filepath, "r") as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    pattern = lines[0]
    text = lines[1]
    d = int(lines[2])

    count = approximate_pattern_count(text, pattern, d)

    if output_filepath:
        with open(output_filepath, "w") as out_file:
            out_file.write(str(count))
        print(f"Results successfully written to {output_filepath}")

    return count


def neighbors(pattern: str, d: int) -> set:
    """Generates all k-mers within Hamming distance d from Pattern (d-neighborhood)."""
    if d == 0:
        return {pattern}
    if len(pattern) == 1:
        return {"A", "C", "G", "T"}

    neighborhood = set()
    suffix_neighbors = neighbors(pattern[1:], d)

    for suffix in suffix_neighbors:
        if hamming_distance(pattern[1:], suffix) < d:
            for nuc in ["A", "C", "G", "T"]:
                neighborhood.add(nuc + suffix)
        else:
            neighborhood.add(pattern[0] + suffix)

    return neighborhood


def frequent_words_with_mismatches(text: str, k: int, d: int) -> str:
    """Finds all most frequent k-mers (with up to d mismatches) using neighborhood expansion."""
    freq_map = {}
    n = len(text)

    # Populate frequency table for every k-mer neighbor
    for i in range(n - k + 1):
        pattern = text[i : i + k]
        neighborhood = neighbors(pattern, d)

        for neighbor in neighborhood:
            if neighbor not in freq_map:
                freq_map[neighbor] = 1
            else:
                freq_map[neighbor] += 1

    # Find maximum count and collect top k-mers
    max_count = maxMap(freq_map)
    frequent_patterns = [
        pattern for pattern, count in freq_map.items() if count == max_count
    ]

    return " ".join(frequent_patterns)


def frequent_words_with_mismatches_from_file(
    input_filepath: str, output_filepath: str = None
) -> str:
    """Reads Text (line 1) and k, d (line 2) from an input file,

    calculates frequent words with mismatches, and optionally writes to file.
    """
    with open(input_filepath, "r") as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    text = lines[0]
    k, d = map(int, lines[1].split())

    result = frequent_words_with_mismatches(text, k, d)

    if output_filepath:
        with open(output_filepath, "w") as out_file:
            out_file.write(result)
        print(f"Results successfully written to {output_filepath}")

    return result


# ==========================================
# Execution Example
# ==========================================
if __name__ == "__main__":
    # 1. Inline Test for Approximate Pattern Count
    sample_text = "AACAAGCTGATAAACATTTAAAGAG"
    sample_pattern = "AAAAA"
    sample_d = 2

    print(
        f"Count2('{sample_text}', '{sample_pattern}') =",
        approximate_pattern_count(sample_text, sample_pattern, sample_d),
    )

    print("\n" + "=" * 50 + "\n")

    # 2. File Execution for Approximate Pattern Count
    filepath_app_count = "dataset_30278_6.txt"
    try:
        count_result = approximate_pattern_count_from_file(filepath_app_count)
        print(
            f"Approximate Pattern Count from '{filepath_app_count}':",
            count_result,
        )
    except FileNotFoundError:
        print(f"File '{filepath_app_count}' not found.")

    print("\n" + "=" * 50 + "\n")

    # 3. Inline Test for Frequent Words with Mismatches
    test_text = (
        "ACGTTGCATGTCACGCTTTATCACGGGACACCCGGCGACACCCGGCGACACCCGGCGACACCCGGCGACGACACCCGGCG"
    )
    k = 4
    d = 1

    result = frequent_words_with_mismatches(test_text, k, d)
    print(
        f"Most frequent {k}-mers with at most {d} mismatch(es) (Neighbor-based):"
    )
    print(result)

    print("\n" + "=" * 50 + "\n")

    # 4. File Execution for Frequent Words with Mismatches
    filepath_fwm = "dataset_30278_9.txt"
    try:
        file_result = frequent_words_with_mismatches_from_file(filepath_fwm)
        print(f"Frequent Words with Mismatches from '{filepath_fwm}':")
        print(file_result)
    except FileNotFoundError:
        print(f"File '{filepath_fwm}' not found.")