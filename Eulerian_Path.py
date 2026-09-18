import sys
from collections import defaultdict
from typing import Iterable


def eulerian_path(g: dict[int, list[int]]) -> Iterable[int]:
    """Constructs an Eulerian path in a directed graph."""
    if not g:
        return []

    # Calculate in-degrees and out-degrees to find start and end nodes
    in_degree = defaultdict(int)
    out_degree = defaultdict(int)

    for u, neighbors in g.items():
        out_degree[u] += len(neighbors)
        for v in neighbors:
            in_degree[v] += 1

    start_node = None
    end_node = None

    # Find nodes where degrees do not balance
    all_nodes = set(g.keys()).union(in_degree.keys())
    for node in all_nodes:
        out_d = out_degree[node]
        in_d = in_degree[node]
        if out_d - in_d == 1:
            start_node = node
        elif in_d - out_d == 1:
            end_node = node

    # If no unbalanced nodes, it's an Eulerian cycle; pick any node with outgoing edges
    if start_node is None:
        start_node = next(iter(g.keys()))

    # Make a copy of graph and add a dummy edge from end_node to start_node if path exists
    adj_list = {u: list(neighbors)[::-1] for u, neighbors in g.items()}
    if end_node is not None and start_node is not None:
        if end_node not in adj_list:
            adj_list[end_node] = []
        adj_list[end_node].append(start_node)

    # Run Hierholzer's Algorithm starting from start_node
    if start_node is None:
        return []
    stack: list[int] = [start_node]
    cycle = []

    while stack:
        curr = stack[-1]
        if adj_list.get(curr):
            next_node = adj_list[curr].pop()
            stack.append(next_node)
        else:
            cycle.append(stack.pop())

    path = cycle[::-1]

    # Remove the added dummy edge (end_node -> start_node) from the resulting cycle
    if end_node is not None:
        for i in range(len(path) - 1):
            if path[i] == end_node and path[i + 1] == start_node:
                path = path[i + 1 :] + path[1 : i + 1]
                break

    return path


def parse_adjacency_list(raw_data_str: str) -> dict[int, list[int]]:
    """Parses colon-formatted string input into an adjacency dictionary."""
    graph = defaultdict(list)
    lines = raw_data_str.strip().splitlines()

    for line in lines:
        line = line.strip()
        if not line:
            continue
        node_str, neighbors_str = line.split(":")
        u = int(node_str.strip())
        neighbors = list(map(int, neighbors_str.strip().split()))
        graph[u] = neighbors

    return dict(graph)


def solve_from_string(raw_data_str: str) -> str:
    """Parses a multi-line graph string, finds the Eulerian path, and returns space-separated nodes."""
    graph = parse_adjacency_list(raw_data_str)
    path = eulerian_path(graph)
    output_str = " ".join(map(str, path))

    print("\nEulerian Path Output:")
    print(output_str)
    return output_str


def solve_from_file(
    input_filename="input.txt", output_filename="output_eulerian_path.txt"
) -> Iterable[int] | None:
    """Reads an adjacency list from an input file, computes the Eulerian path, prints it, and writes the output to a file."""
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        graph = parse_adjacency_list(content)
        path = eulerian_path(graph)
        output_str = " ".join(map(str, path))

        print(
            f"\n[+] Successfully processed Eulerian Path ({len(graph)} nodes):"
        )
        print("-" * 40)
        print(output_str)
        print("-" * 40)

        with open(output_filename, "w") as f:
            f.write(output_str + "\n")

        print(f"[+] Output written to '{output_filename}'")
        return path

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


def solve_from_cmd():
    """Reads input directly from standard input (stdin/cmd pipe)."""
    content = sys.stdin.read().strip()
    if not content:
        return
    graph = parse_adjacency_list(content)
    path = eulerian_path(graph)
    print(" ".join(map(str, path)))


if __name__ == "__main__":
    if not sys.stdin.isatty():
        solve_from_cmd()
    elif len(sys.argv) > 1:
        solve_from_file(input_filename=sys.argv[1])
    else:
        print("Eulerian Path Solver")
        print("=" * 35)

        sample_input = """0: 2
1: 3
2: 1
3: 0 4
6: 3 7
7: 8
8: 9
9: 6"""
        print(f"\nRunning Sample Input 1:\n{sample_input}\n")
        solve_from_string(sample_input)