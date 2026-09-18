import sys
from collections import defaultdict

# Increase Python's maximum recursion limit for large Eulerian paths
sys.setrecursionlimit(200000)


def build_paired_de_bruijn_graph(paired_reads: list[tuple[str, str]]):
    """Constructs a paired de Bruijn graph along with in-degrees and out-degrees."""
    graph = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for kmer1, kmer2 in paired_reads:
        prefix = (kmer1[:-1], kmer2[:-1])
        suffix = (kmer1[1:], kmer2[1:])
        graph[prefix].append(suffix)
        out_degree[prefix] += 1
        in_degree[suffix] += 1

    return dict(graph), in_degree, out_degree


def string_spelled_by_gapped_patterns(
    gapped_patterns: list[tuple[str, str]], k: int, d: int
) -> str | None:
    """
    Reconstructs a genome string from an ordered list of (k-1)-mer pairs.
    Assembles FirstPatterns and SecondPatterns and checks if PrefixString and SuffixString
    overlap consistently at offset (k + d).
    """
    if not gapped_patterns:
        return None

    first_patterns = [pair[0] for pair in gapped_patterns]
    second_patterns = [pair[1] for pair in gapped_patterns]

    prefix_string = first_patterns[0] + "".join(kmer[-1] for kmer in first_patterns[1:])
    suffix_string = second_patterns[0] + "".join(kmer[-1] for kmer in second_patterns[1:])

    overlap_len = len(prefix_string) - (k + d)
    if prefix_string[k + d :] == suffix_string[:overlap_len]:
        return prefix_string + suffix_string[overlap_len:]

    return None


def find_valid_eulerian_path(
    graph: dict[tuple[str, str], list[tuple[str, str]]],
    in_degree: dict[tuple[str, str], int],
    out_degree: dict[tuple[str, str], int],
    k: int,
    d: int,
    num_edges: int,
) -> str | None:
    """Finds an Eulerian path in the paired de Bruijn graph spelling a valid genome string."""
    all_nodes = set(graph.keys()).union(in_degree.keys())

    start_nodes = []
    for node in all_nodes:
        if out_degree[node] - in_degree[node] == 1:
            start_nodes = [node]
            break

    if not start_nodes:
        start_nodes = list(graph.keys())

    for start_node in start_nodes:
        adj_list = {u: list(neighbors) for u, neighbors in graph.items()}
        path = [start_node]

        def backtrack(curr_node):
            if len(path) == num_edges + 1:
                spelled_string = string_spelled_by_gapped_patterns(path, k - 1, d + 1)
                if spelled_string is not None:
                    return spelled_string
                return None

            neighbors = adj_list.get(curr_node, [])
            for i in range(len(neighbors)):
                next_node = neighbors.pop(i)
                path.append(next_node)

                res = backtrack(next_node)
                if res is not None:
                    return res

                path.pop()
                neighbors.insert(i, next_node)

            return None

        result = backtrack(start_node)
        if result is not None:
            return result

    return None


def string_reconstruction_read_pairs(arg1, arg2, arg3=None) -> str:
    """
    Solves the String Reconstruction from Read-Pairs Problem.
    Handles flexible argument ordering: (paired_reads, k, d) or (k, d, paired_reads).
    """
    if isinstance(arg1, (list, tuple)):
        paired_reads, k, d = arg1, arg2, arg3
    elif isinstance(arg3, (list, tuple)):
        k, d, paired_reads = arg1, arg2, arg3
    else:
        paired_reads, k, d = arg1, arg2, arg3

    if paired_reads and isinstance(paired_reads[0], str):
        paired_reads = [tuple(read.split("|")) for read in paired_reads]

    if not isinstance(k, int) or not isinstance(d, int):
        raise TypeError("k and d must be integers")

    num_edges = len(paired_reads)
    graph, in_degree, out_degree = build_paired_de_bruijn_graph(paired_reads) # pyright: ignore[reportArgumentType]

    result = find_valid_eulerian_path(
        graph, in_degree, out_degree, k, d, num_edges
    )
    return result if result is not None else ""


# --- INPUT WRAPPER FUNCTIONS --- #


def solve_from_string(raw_data_str: str) -> str:
    """Parses a raw whitespace-separated input string and solves the problem."""
    tokens = raw_data_str.strip().split()
    if not tokens:
        return ""
    k = int(tokens[0])
    d = int(tokens[1])
    paired_reads = [tuple(token.split("|")) for token in tokens[2:]]
    return string_reconstruction_read_pairs(k, d, paired_reads)


def solve_from_file(filename: str) -> str:
    """Reads input data from a specified file path and returns the solution."""
    with open(filename, "r") as f:
        content = f.read()
    return solve_from_string(content)


def solve_from_cmd() -> str:
    """Reads input data from standard command line input (sys.stdin)."""
    content = sys.stdin.read()
    return solve_from_string(content)


if __name__ == "__main__":
    # If piped input or redirected stdin is provided
    if not sys.stdin.isatty():
        print(solve_from_cmd())
    # If a filename argument is passed: python script.py dataset.txt
    elif len(sys.argv) > 1:
        print(solve_from_file(sys.argv[1]))
    # Fallback default sample test
    else:
        sample_input = "4 2 GTTT|ATTT TTTA|TTTG TTAC|TTGT TACG|TGTA ACGT|GTAT CGTT|TATT"
        print("Reconstructed String:")
        print(solve_from_string(sample_input))