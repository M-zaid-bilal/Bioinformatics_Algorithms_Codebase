import sys
from collections import defaultdict
from typing import List, Dict,Iterable


def eulerian_cycle(graph: dict[int, list[int]]) -> list[int]:
    """
    Finds an Eulerian cycle in a graph using Hierholzer's Algorithm in O(E) time.
    """
    if not graph:
        return []

    # Deep copy graph to avoid mutating input & reverse for O(1) popping
    adj_list = {u: neighbors[::-1] for u, neighbors in graph.items()}
    
    # Start at the first available node key
    start_node = next(iter(adj_list.keys()))
    
    stack = [start_node]
    cycle = []

    while stack:
        curr = stack[-1]
        if adj_list.get(curr):
            next_node = adj_list[curr].pop()
            stack.append(next_node)
        else:
            cycle.append(stack.pop())

    return cycle[::-1]


def parse_adjacency_list(raw_data_str: str) -> dict[int, list[int]]:
    """Parses colon-formatted string input into an adjacency dictionary."""
    graph = defaultdict(list)
    lines = raw_data_str.strip().splitlines()
    
    for line in lines:
        line = line.strip()
        if not line:
            continue
        node_str, neighbors_str = line.split(':')
        u = int(node_str.strip())
        neighbors = list(map(int, neighbors_str.strip().split()))
        graph[u] = neighbors
        
    return dict(graph)


def solve_from_string(raw_data_str: str) -> str:
    """Parses a multi-line graph string, finds the Eulerian cycle, and returns space-separated nodes."""
    graph = parse_adjacency_list(raw_data_str)
    cycle = eulerian_cycle(graph)
    output_str = " ".join(map(str, cycle))
    
    print("\nEulerian Cycle Output:")
    print(output_str)
    return output_str


def solve_from_file(input_filename="input_6_ec.txt", output_filename="output_eulerian_cycle.txt") -> list[int] | None:
    """
    Reads an adjacency list from an input file, computes the Eulerian cycle,
    prints it, and writes the output to a file.
    """
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        graph = parse_adjacency_list(content)
        cycle = eulerian_cycle(graph)
        output_str = " ".join(map(str, cycle))

        print(f"\n[+] Successfully processed Eulerian Cycle ({len(graph)} nodes):")
        print("-" * 40)
        print(output_str)
        print("-" * 40)

        with open(output_filename, "w") as f:
            f.write(output_str + "\n")

        print(f"[+] Output written to '{output_filename}'")
        return cycle

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


def solve_from_cmd():
    """Reads input directly from standard input (stdin/cmd pipe)."""
    content = sys.stdin.read().strip()
    if not content:
        return
    graph = parse_adjacency_list(content)
    cycle = eulerian_cycle(graph)
    print(" ".join(map(str, cycle)))


if __name__ == "__main__":
    # If piped data or file args exist, handle accordingly; otherwise run sample demo
    if not sys.stdin.isatty():
        solve_from_cmd()
    elif len(sys.argv) > 1:
        solve_from_file(input_filename=sys.argv[1])
    else:
        print("Eulerian Cycle Solver")
        print("=" * 35)

        sample_input = """0: 3
1: 0
2: 1 6
3: 2
4: 2
5: 4
6: 5 8
7: 9
8: 7
9: 6"""
        print(f"\nRunning Sample Input:\n{sample_input}\n")
        solve_from_string(sample_input)