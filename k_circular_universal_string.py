import sys
# Import functions directly from your debrujin_graph.py module
from debrujin_graph import build_de_bruijn_graph, find_eulerian_cycle


def k_universal_circular_string(k: int) -> str:
    """
    Generates a k-universal binary circular string of length 2^k
    using imported functions from debrujin_graph.
    """
    if k == 1:
        return "01"

    graph = build_de_bruijn_graph(k)
    start_node = "0" * (k - 1)

    path_nodes = find_eulerian_cycle(graph, start_node)

    # Extract the last character of each node in the cycle traversal
    circular_string = "".join(node[-1] for node in path_nodes[1:])
    return circular_string


def parse_input_string(raw_data_str: str) -> int:
    """Parses whitespace-separated string into integer k."""
    tokens = raw_data_str.strip().split()
    if not tokens:
        return 0
    return int(tokens[0])


def solve_from_string(raw_data_str: str) -> str:
    """Parses input string, runs k-universal circular string generator, and prints result."""
    k = parse_input_string(raw_data_str)
    result = k_universal_circular_string(k)

    print(f"\n{k}-Universal Circular String:")
    print(result)
    return result


def solve_from_file(
    input_filename="input.txt", output_filename="output_k_universal.txt"
) -> str | None:
    """Reads k from file, computes the k-universal circular string, prints it, and writes output to file."""
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        k = parse_input_string(content)
        result = k_universal_circular_string(k)

        print(f"\n[+] Successfully generated {k}-Universal Circular String (length {len(result)}):")
        print("-" * 40)
        print(result)
        print("-" * 40)

        with open(output_filename, "w") as f:
            f.write(result + "\n")

        print(f"[+] Output written to '{output_filename}'")
        return result

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


def solve_from_cmd():
    """Reads k directly from standard input (stdin/cmd pipe)."""
    content = sys.stdin.read().strip()
    if not content:
        return
    k = parse_input_string(content)
    result = k_universal_circular_string(k)
    print(result)


if __name__ == "__main__":
    if not sys.stdin.isatty():
        solve_from_cmd()
    elif len(sys.argv) > 1:
        solve_from_file(input_filename=sys.argv[1])
    else:
        print("k-Universal Circular String Solver")
        print("=" * 35)

        sample_input = "3"
        print(f"\nRunning Sample Input:\n{sample_input}\n")
        solve_from_string(sample_input)