import sys
from collections import defaultdict


def build_de_bruijn_graph(kmers: list[str]):
    """Constructs a De Bruijn graph from k-mers along with in/out degree counts."""
    graph = defaultdict(list)
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)
    all_nodes = set()

    for kmer in kmers:
        prefix = kmer[:-1]
        suffix = kmer[1:]
        graph[prefix].append(suffix)
        out_degree[prefix] += 1
        in_degree[suffix] += 1
        all_nodes.add(prefix)
        all_nodes.add(suffix)

    return graph, in_degree, out_degree, all_nodes


def get_maximal_non_branching_paths(kmers: list[str]) -> list[list[str]]:
    """Finds all maximal non-branching paths (contigs) in the De Bruijn graph."""
    graph, in_degree, out_degree, all_nodes = build_de_bruijn_graph(kmers)
    paths = []
    visited_edges = set()

    def is_one_in_one_out(v):
        return in_degree[v] == 1 and out_degree[v] == 1

    # 1. Find maximal non-branching paths starting from non-1-in-1-out nodes
    for v in all_nodes:
        if not is_one_in_one_out(v):
            if out_degree[v] > 0:
                for w in graph[v]:
                    path = [v, w]
                    visited_edges.add((v, w))
                    curr = w
                    while is_one_in_one_out(curr):
                        next_node = graph[curr][0]
                        path.append(next_node)
                        visited_edges.add((curr, next_node))
                        curr = next_node
                    paths.append(path)

    # 2. Find isolated cycles consisting entirely of 1-in-1-out nodes
    for v in all_nodes:
        if is_one_in_one_out(v):
            # Check if any outgoing edge from v has not been visited yet
            for w in graph[v]:
                if (v, w) not in visited_edges:
                    cycle = [v, w]
                    visited_edges.add((v, w))
                    curr = w
                    while curr != v:
                        next_node = graph[curr][0]
                        cycle.append(next_node)
                        visited_edges.add((curr, next_node))
                        curr = next_node
                    paths.append(cycle)

    return paths


def generate_contigs(kmers: list[str]) -> list[str]:
    """Generates contiguous sequence strings from maximal non-branching paths."""
    paths = get_maximal_non_branching_paths(kmers)
    contigs = []

    for path in paths:
        # Reconstruct sequence: first node + remaining character of subsequent nodes
        contig = path[0] + "".join(node[-1] for node in path[1:])
        contigs.append(contig)

    return contigs


# --- INPUT WRAPPER FUNCTIONS --- #


def solve_from_string(raw_data_str: str) -> str:
    """Parses whitespace-separated k-mers and returns space-separated contigs."""
    kmers = raw_data_str.strip().split()
    if not kmers:
        return ""
    contigs = generate_contigs(kmers)
    result = " ".join(sorted(contigs))
    return result


def solve_from_file(filename: str) -> str:
    """Reads k-mers from a file and returns contigs."""
    with open(filename, "r") as f:
        content = f.read()
    return solve_from_string(content)


def solve_from_cmd() -> str:
    """Reads k-mers from standard input."""
    content = sys.stdin.read()
    return solve_from_string(content)


if __name__ == "__main__":
    if not sys.stdin.isatty():
        print(solve_from_cmd())
    elif len(sys.argv) > 1:
        print(solve_from_file(sys.argv[1]))
    else:
        sample_input = "ATG ATG TGT TGG CAT GGA GAT AGA"
        print("Generated Contigs:")
        print(solve_from_string(sample_input))