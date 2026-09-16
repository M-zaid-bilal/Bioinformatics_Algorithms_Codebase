def genome_path_to_string(patterns: list) -> str:
    """
    Reconstructs a genome string from a path of k-mers.
    
    Parameters:
    - patterns (list of str): Ordered sequence of k-mers where the last (k-1) 
                              symbols of Pattern[i] equal the first (k-1) 
                              symbols of Pattern[i+1].
    
    Returns:
    - str: The reconstructed genome sequence of length k + n - 1.
    """
    if not patterns:
        return ""
    
    # Start with the full first k-mer
    reconstructed = patterns[0]
    
    # Append the last character of each subsequent k-mer in the path
    for pattern in patterns[1:]:
        reconstructed += pattern[-1]
        
    return reconstructed



def solve_from_file(input_filename="dataset_genome_path.txt", output_filename="output_genome_path.txt"):
    """
    Reads k-mer patterns from an input dataset file, reconstructs the genome string,
    and writes the output to a file. Handles both space-separated and line-separated inputs.
    """
    try:
        with open(input_filename, "r") as f:
            content = f.read().strip()

        if not content:
            print(f"Error: File '{input_filename}' is empty.")
            return None

        # Split on any whitespace to handle both single-line and multi-line datasets
        patterns = content.split()

        result = genome_path_to_string(patterns)

        with open(output_filename, "w") as f:
            f.write(result + "\n")

        print(f"\n[+] Successfully processed {len(patterns)} k-mers.")
        print(f"[+] Reconstructed string length: {len(result)}")
        print(f"[+] Reconstructed String:\n{result}")
        print(f"[+] Output written to '{output_filename}'")
        return result

    except FileNotFoundError:
        print(f"Error: Input file '{input_filename}' not found.")
        return None



def solve_from_string(raw_data_str: str) -> str:
    """
    Parses a raw multi-line or space-separated string (e.g., copied from Stepik/Cogniterra)
    and prints/returns the reconstructed genome sequence.
    """
    patterns = raw_data_str.strip().split()
    result = genome_path_to_string(patterns)
    print("\nReconstructed Genome Path String:")
    print(result)
    return result



if __name__ == "__main__":
    # Example usage / interactive runner
    print("Genome Path to String Solver")
    print("-" * 30)
    
    choice = input("Enter '1' for File Mode, or '2' to paste Raw Input: ").strip()
    
    if choice == '1':
        filename = input("Enter input filename (or press Enter for 'input_1_gp.txt'): ").strip()
        if not filename:
            filename = "input_1_gp.txt"
        solve_from_file(filename)
    elif choice == '2':
        print("\nPaste your space-separated or multi-line k-mers below:")
        raw_input = input().strip()
        solve_from_string(raw_input)
    else:
        # Default run with Sample 1
        sample_input = "ACCGA CCGAA CGAAG GAAGC AAGCT"
        print(f"\nRunning default Sample Input: {sample_input}")
        solve_from_string(sample_input)