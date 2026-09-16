import sys
def de_bruijn_from_string(text: str, k: int) -> dict[str, list[str]]:
    """
    Constructs the De Bruijn graph DeBruijn_k(Text) from a given string.
    
    Nodes in DeBruijn_k(Text) correspond to all unique (k-1)-mers present in Text.
    Directed edges correspond to all overlapping k-mers in Text, pointing from
    the (k-1)-mer prefix to the (k-1)-mer suffix.
    """
    adj_matrix: dict[str, list[str]] = {}

    # Iterate through all overlapping k-mers in Text
    for i in range(len(text) - k + 1):
        kmer = text[i : i + k]
        prefix = kmer[:-1]  # Left (k-1)-mer
        suffix = kmer[1:]   # Right (k-1)-mer

        if prefix not in adj_matrix:
            adj_matrix[prefix] = []
        adj_matrix[prefix].append(suffix)

    # Sort target nodes lexicographically for each source node
    for prefix in adj_matrix:
        adj_matrix[prefix].sort()

    return adj_matrix


def format_adjacency_list(adj_dict: dict[str, list[str]], style: str = "colon") -> list[str]:
    """
    Formats the adjacency dictionary into standard representation lines.
    
    Styles supported:
      - 'colon' (default):  'GT: TA TG'
      - 'arrow':          'GT -> TA, TG'
      - 'arrow_space':    'GT -> TA,TG'
    """
    # Sort source nodes lexicographically
    sorted_nodes = sorted(adj_dict.keys())
    lines = []

    for u in sorted_nodes:
        targets = adj_dict[u]
        if style == "colon":
            target_str = " ".join(targets)
            lines.append(f"{u}: {target_str}")
        elif style == "arrow":
            target_str = ", ".join(targets)
            lines.append(f"{u} -> {target_str}")
        elif style == "arrow_space":
            target_str = ",".join(targets)
            lines.append(f"{u} -> {target_str}")

    return lines


def solve_from_string(raw_data_str: str, style: str = "colon") -> str:
    """
    Parses raw string input containing k on line 1 and Text on line 2,
    then outputs the De Bruijn Graph adjacency list.
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    if len(lines) < 2:
        print("Error: Input must contain k on line 1 and Text on line 2.")
        return ""

    k = int(lines[0])
    text = lines[1]

    adj_dict = de_bruijn_from_string(text, k)
    formatted_lines = format_adjacency_list(adj_dict, style=style)
    output_str = "\n".join(formatted_lines)

    print(output_str)
    return output_str


def solve_from_file(input_filename="input_1_dbs.txt", output_filename="output_debruijn_string.txt", style="colon"):
    """
    Reads k and Text from an input dataset file, computes DeBruijn_k(Text),
    and writes the resulting adjacency list to an output file.
    """
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        lines = [line.strip() for line in content.split("\n") if line.strip()]
        k = int(lines[0])
        text = lines[1]

        adj_dict = de_bruijn_from_string(text, k)
        formatted_lines = format_adjacency_list(adj_dict, style=style)
        output_str = "\n".join(formatted_lines)

        with open(output_filename, "w") as f:
            f.write(output_str + "\n")

        print(f"\n[+] De Bruijn Graph DeBruijn_{k}(Text) generated successfully ({len(adj_dict)} nodes):")
        print("-" * 40)
        print(output_str)
        print("-" * 40)
        print(f"[+] Output written to '{output_filename}'")
        return output_str

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


if __name__ == "__main__":
    # Test Sample 1
    sample_input = """3
ACGTGTATA"""

    print("=== Testing Sample 1 ===")
    solve_from_string(sample_input, style="colon")