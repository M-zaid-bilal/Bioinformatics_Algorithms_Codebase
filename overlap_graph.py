def overlap_graph(patterns: list[str]) -> dict[str, list[str]]:
    """Forms the overlap graph of a collection of patterns."""
    adj_list: dict[str, list[str]] = {}
    
    sorted_patterns = sorted(patterns)
    
    for i, u in enumerate(sorted_patterns):
        suffix = u[1:]
        matches = []
        for j, v in enumerate(sorted_patterns):
            # Allow self-loop overlaps (e.g., TTT -> TTT) where Suffix(u) == Prefix(v)
            if suffix == v[:-1]:
                matches.append(v)
        if matches:
            if u not in adj_list:
                adj_list[u] = []
            adj_list[u].extend(matches)
            
    for u in adj_list:
        adj_list[u] = sorted(list(dict.fromkeys(adj_list[u])))
            
    return adj_list


def format_adjacency_list(adj_dict: dict, style="colon") -> list:
    """
    Formats the adjacency dictionary into standard string representations:
    - 'colon':  'AAG: AGA' or 'TCT: CTA CTC'
    - 'arrow':  'AAG -> AGA'
    - 'comma':  'TCT -> CTA, CTC'
    """
    lines = []
    for u, targets in adj_dict.items():
        if style == "colon":
            target_str = " ".join(targets)
            lines.append(f"{u}: {target_str}")
        elif style == "comma":
            target_str = ", ".join(targets)
            lines.append(f"{u} -> {target_str}")
        else:  # arrow per line
            for v in targets:
                lines.append(f"{u} -> {v}")
    return lines


def solve_from_file(input_filename="dataset_overlap_graph.txt", output_filename="output_overlap_graph.txt", style="colon"):
    """
    Reads k-mer patterns from an input file, constructs the overlap graph,
    prints the result to the console, and writes the output to a file.
    """
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: Input file '{input_filename}' is empty.")
            return None

        # Split on any whitespace to handle both line-separated and space-separated inputs
        patterns = content.split()

        adj_dict = overlap_graph(patterns)
        formatted_lines = format_adjacency_list(adj_dict, style=style)
        output_str = "\n".join(formatted_lines)

        # Print directly to console
        print(f"\n[+] Successfully constructed Overlap Graph ({len(patterns)} k-mers, {len(formatted_lines)} nodes):")
        print("-" * 40)
        print(output_str)
        print("-" * 40)

        # Save output to file
        with open(output_filename, "w") as f:
            f.write(output_str + "\n")

        print(f"[+] Output written to '{output_filename}'")
        return adj_dict

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None


def solve_from_string(raw_data_str: str, style="colon") -> str:
    """
    Parses a raw space-separated or multi-line string dataset (e.g., copied from Stepik/Cogniterra),
    prints the constructed overlap graph, and returns the formatted result.
    """
    patterns = raw_data_str.strip().split()
    adj_dict = overlap_graph(patterns)
    formatted_lines = format_adjacency_list(adj_dict, style=style)
    output_str = "\n".join(formatted_lines)

    print("\nOverlap Graph Adjacency List:")
    print(output_str)
    return output_str



if __name__ == "__main__":
    print("Overlap Graph Construction Solver")
    print("=" * 35)

    sample_input = "AAG AGA ATT CTA CTC GAT TAC TCT TCT TTC"
    print(f"\nRunning Sample Input 1:\n{sample_input}\n")
    solve_from_string(sample_input, style="colon")