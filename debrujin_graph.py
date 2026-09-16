import itertools


def build_de_bruijn_graph(k: int) -> dict[str, list[str]]:
    """
    Constructs a De Bruijn graph dB_k for binary k-mers.
    Nodes are binary (k-1)-mers, and directed edges correspond to k-mers.
    """
    nodes = ["".join(p) for p in itertools.product("01", repeat=k - 1)]
    graph: dict[str, list[str]] = {node: [] for node in nodes}

    for node in nodes:
        # Every (k-1)-mer u connects to u[1:] + '0' and u[1:] + '1'
        suffix = node[1:]
        graph[node].append(suffix + "0")
        graph[node].append(suffix + "1")

    return graph


def find_eulerian_cycle(graph: dict[str, list[str]], start_node: str) -> list[str]:
    """
    Finds an Eulerian cycle in a directed Eulerian graph using Hierholzer's Algorithm.
    """
    # Create a mutable copy of adjacency lists
    adj = {u: list(v_list) for u, v_list in graph.items()}
    
    stack = [start_node]
    cycle = []

    while stack:
        u = stack[-1]
        if adj[u]:
            v = adj[u].pop()
            stack.append(v)
        else:
            cycle.append(stack.pop())

    # Cycle is reconstructed in reverse order
    return cycle[::-1]


def k_universal_string(k: int, circular: bool = True) -> str:
    """
    Generates a k-universal binary string of length 2^k (circular)
    or length 2^k + k - 1 (linear/non-circular).
    """
    if k == 1:
        return "01" if circular else "010"

    graph = build_de_bruijn_graph(k)
    start_node = "0" * (k - 1)
    
    path_nodes = find_eulerian_cycle(graph, start_node)

    # Extract the last character of each node along the Eulerian path traversal
    circular_string = "".join(node[-1] for node in path_nodes[1:])

    if circular:
        return circular_string
    else:
        # Extend with the first (k-1) characters to form a linear string
        return circular_string + circular_string[: k - 1]


def verify_k_universal(universal_str: str, k: int, circular: bool = True) -> bool:
    """
    Validates whether a string contains all 2^k unique binary k-mers.
    """
    expected_kmers = {"".join(p) for p in itertools.product("01", repeat=k)}
    found_kmers = set()

    if circular:
        extended_str = universal_str + universal_str[: k - 1]
    else:
        extended_str = universal_str

    for i in range(len(extended_str) - k + 1):
        found_kmers.add(extended_str[i : i + k])

    return found_kmers == expected_kmers


if __name__ == "__main__":
    k = 4
    
    circular_solution = k_universal_string(k, circular=True)
    linear_solution = k_universal_string(k, circular=False)

    print(f"=== {k}-Universal String Solution ===")
    print(f"Circular String  (length {len(circular_solution)}):  {circular_solution}")
    print(f"Linear String    (length {len(linear_solution)}):  {linear_solution}")
    
    is_valid = verify_k_universal(circular_solution, k, circular=True)
    print(f"\nVerification status: {'PASS (Contains all 16 4-mers exactly once)' if is_valid else 'FAIL'}")