import sys
from collections import defaultdict
from typing import Iterable


def de_bruijn_from_patterns(patterns: list[str]) -> dict[str, list[str]]:
    """Constructs a de Bruijn graph from a collection of k-mers."""
    graph = defaultdict(list)
    for pattern in patterns:
        prefix = pattern[:-1]
        suffix = pattern[1:]
        graph[prefix].append(suffix)
    return dict(graph)


def eulerian_path(g: dict[str, list[str]]) -> list[str]:
    """Constructs an Eulerian path in a directed graph using Hierholzer's Algorithm."""
    if not g:
        return []

    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for u, neighbors in g.items():
        out_degree[u] += len(neighbors)
        for v in neighbors:
            in_degree[v] += 1

    start_node = None
    end_node = None

    all_nodes = set(g.keys()).union(in_degree.keys())
    for node in all_nodes:
        out_d = out_degree[node]
        in_d = in_degree[node]
        if out_d - in_d == 1:
            start_node = node
        elif in_d - out_d == 1:
            end_node = node

    if start_node is None:
        start_node = next(iter(g.keys()))
    if start_node is None:
        return []

    adj_list = {u: list(neighbors)[::-1] for u, neighbors in g.items()}
    if end_node is not None and start_node is not None:
        if end_node not in adj_list:
            adj_list[end_node] = []
        adj_list[end_node].append(start_node)

    stack = [start_node]
    cycle = []

    while stack:
        curr = stack[-1]
        if curr is not None and adj_list.get(curr):
            next_node = adj_list[curr].pop()
            stack.append(next_node)
        else:
            cycle.append(stack.pop())

    path = cycle[::-1]

    if end_node is not None:
        for i in range(len(path) - 1):
            if path[i] == end_node and path[i + 1] == start_node:
                path = path[i + 1 :] + path[1 : i + 1]
                break

    return path


def path_to_genome(path: list[str]) -> str:
    """Reconstructs a genome string from an ordered list of (k-1)-mers."""
    if not path:
        return ""
    genome = path[0]
    for kmer in path[1:]:
        genome += kmer[-1]
    return genome


def string_reconstruction(k: int, patterns: list[str]) -> str:
    """Solves the String Reconstruction Problem via de Bruijn Graph & Eulerian Path."""
    db_graph = de_bruijn_from_patterns(patterns)
    path = eulerian_path(db_graph)
    return path_to_genome(path)


def parse_input_string(raw_data_str: str) -> tuple[int, list[str]]:
    """Parses whitespace-separated string into integer k and list of patterns."""
    tokens = raw_data_str.strip().split()
    if not tokens:
        return 0, []
    k = int(tokens[0])
    patterns = tokens[1:]
    return k, patterns


def solve_from_string(raw_data_str: str) -> str:
    """Parses input string, runs string reconstruction, and prints output."""
    k, patterns = parse_input_string(raw_data_str)
    text = string_reconstruction(k, patterns)

    print("\nReconstructed Genome String:")
    print(text)
    return text


def solve_from_file(
    input_filename="input.txt", output_filename="output_string_reconstruction.txt"
) -> str | None:
    """Reads dataset from a file, reconstructs the genome, prints it, and writes output to a file."""
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        k, patterns = parse_input_string(content)
        text = string_reconstruction(k, patterns)

        print(
            f"\n[+] Successfully reconstructed string of length {len(text)} from {len(patterns)} {k}-mers:"
        )
        print("-" * 40)
        print(text)
        print("-" * 40)

        with open(output_filename, "w") as f:
            f.write(text + "\n")

        print(f"[+] Output written to '{output_filename}'")
        return text

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


def solve_from_cmd():
    """Reads input directly from standard input (stdin/cmd pipe)."""
    content = sys.stdin.read().strip()
    if not content:
        return
    k, patterns = parse_input_string(content)
    text = string_reconstruction(k, patterns)
    print(text)


if __name__ == "__main__":
    if not sys.stdin.isatty():
        solve_from_cmd()
    elif len(sys.argv) > 1:
        solve_from_file(input_filename=sys.argv[1])
    else:
        print("String Reconstruction Solver")
        print("=" * 35)

        sample_input = "4 CTAC CTCC TCCT ACTC CCTC CCTA TACT"
        print(f"\nRunning Sample Input:\n{sample_input}\n")
        solve_from_string(sample_input)