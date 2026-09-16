import sys
from collections import defaultdict


def de_bruijn_from_kmers(patterns: list[str]) -> dict[str, list[str]]:
    """
    Constructs the De Bruijn graph DeBruijn(Patterns) from a collection of k-mers.
    
    Nodes in DeBruijn(Patterns) correspond to all unique (k-1)-mers extracted as 
    prefixes and suffixes of the k-mers in Patterns.
    Directed edges correspond to the k-mers themselves, pointing from kmer[:-1] to kmer[1:].
    """
    adj_matrix: dict[str, list[str]] = defaultdict(list)

    for kmer in patterns:
        prefix = kmer[:-1]  # Left (k-1)-mer
        suffix = kmer[1:]   # Right (k-1)-mer
        adj_matrix[prefix].append(suffix)

    for prefix in adj_matrix:
        adj_matrix[prefix].sort()

    return adj_matrix


def format_adjacency_list(adj_dict: dict[str, list[str]], style: str = "colon") -> list[str]:
    """
    Formats the adjacency dictionary into standard representation lines.
    
    Styles supported:
      - 'colon' (default):  'CAG: AGG AGG'
      - 'arrow':          'CAG -> AGG, AGG'
      - 'arrow_space':    'CAG -> AGG,AGG'
    """
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
    Parses space-separated or newline-separated k-mers from a raw string,
    computes DeBruijn(Patterns), and returns the formatted output.
    """
    patterns = raw_data_str.strip().split()
    if not patterns:
        print("Error: Input contains no k-mers.")
        return ""

    adj_dict = de_bruijn_from_kmers(patterns)
    formatted_lines = format_adjacency_list(adj_dict, style=style)
    output_str = "\n".join(formatted_lines)

    print(output_str)
    return output_str


def solve_from_file(input_filename="dataset_debruijn_kmers.txt", output_filename="output_debruijn_kmers.txt", style="colon"):
    """
    Reads k-mer patterns from an input file, constructs DeBruijn(Patterns),
    and writes the resulting adjacency list to an output file.
    """
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        patterns = content.split()
        adj_dict = de_bruijn_from_kmers(patterns)
        formatted_lines = format_adjacency_list(adj_dict, style=style)
        output_str = "\n".join(formatted_lines)

        with open(output_filename, "w") as f:
            f.write(output_str + "\n")

        print(f"\n[+] De Bruijn Graph DeBruijn(Patterns) generated successfully ({len(adj_dict)} source nodes):")
        print("-" * 40)
        print(output_str)
        print("-" * 40)
        print(f"[+] Output written to '{output_filename}'")
        return output_str

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Allows running directly like: python de_bruijn_from_kmers.py my_dataset.txt
        user_filename = sys.argv[1]
        solve_from_file(user_filename)
    else:
        # Interactive prompt if no command-line argument is supplied
        filename_input = input("Enter input dataset filename (press Enter for default 'dataset_debruijn_kmers.txt'): ").strip()
        target_file = filename_input if filename_input else "dataset_debruijn_kmers.txt"
        solve_from_file(target_file)