def hamming_distance(p: str, q: str) -> int:
    """Calculates the Hamming distance between two strings of equal length."""
    distance = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            distance += 1
    return distance


def approximate_pattern_matching(text: str, pattern: str, d: int) -> list[int]:
    """Finds all starting positions of Pattern in Text with at most d mismatches."""
    positions = []
    pattern_len = len(pattern)

    # Slide a window of length `pattern_len` across `text`
    for i in range(len(text) - pattern_len + 1):
        substring = text[i : i + pattern_len]
        if hamming_distance(pattern, substring) <= d:
            positions.append(i)

    return positions


def approximate_pattern_matching_from_file(
    genome_filepath: str, output_filepath: str = None
) -> list[int]:
    """Reads Pattern, Text, and d from a file, then finds matching positions."""
    with open(genome_filepath, "r") as file:
        lines = [line.strip() for line in file.readlines() if line.strip()]

    # Extract inputs according to the file format
    pattern = lines[0]
    text = lines[1]
    d = int(lines[2])

    positions = approximate_pattern_matching(text, pattern, d)

    if output_filepath:
        with open(output_filepath, "w") as out_file:
            out_file.write(" ".join(map(str, positions)))
        print(f"Results successfully written to {output_filepath}")

    return positions


if __name__ == "__main__":
    # Inline Example
    text = "CGCCCGAATCCAGAACGCATTCCCATGTACACACCATACCCCTCCAGCCACCACCCACCACACCCACACACCCACAGCCACCACCCACCACACCCACACACCCACAGCCACCACCCACCACACCCACACACCCAC"
    pattern = "ATTCTGGA"
    d = 3

    positions = approximate_pattern_matching(text, pattern, d)

    print("Pattern:", pattern)
    print("Allowed mismatches (d):", d)
    print("Space separated starting positions:")
    print(" ".join(map(str, positions)))
    print("\n" + "=" * 50 + "\n")

    # File Execution Example
    filepath = "dataset_30278_4.txt"
    file_positions = approximate_pattern_matching_from_file(filepath)
    print(f"Results from file '{filepath}':")
    print(" ".join(map(str, file_positions)))
    string1= "CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA"
    string2= "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"
    print("Hamming Distance : ")
    print(hamming_distance(string1,string2))