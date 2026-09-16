import sys


def kmer_composition(text: str, k: int):
    """
    Generates the k-mer composition of a string Text using a sliding window.
    
    Parameters:
    - text (str): The genomic sequence or text.
    - k (int): Length of each k-mer window.
    
    Returns:
    - list of str: All overlapping k-mers of length k in order of appearance.
    """
    return [text[i:i+k] for i in range(len(text) - k + 1)]



def solve_from_file(input_filename="dataset_string_composition.txt", output_filename="output_string_composition.txt"):
    """
    Reads k and Text from an input dataset file, computes Composition_k(Text),
    and writes each k-mer on a new line to an output file.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        if len(lines) < 2:
            print("Error: Input file must contain 'k' on line 1 and 'Text' on line 2.")
            return None

        k = int(lines[0])
        text = lines[1]

        kmers = kmer_composition(text, k)

        with open(output_filename, "w") as f:
            for kmer in kmers:
                f.write(kmer + "\n")

        print(f"\n[+] Generated {len(kmers)} k-mers of length {k}.")
        print(f"[+] Output written successfully to '{output_filename}'.")
        return kmers

    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None



def solve_from_string(raw_data_str):
    """
    Parses a raw multi-line string (e.g. copied from Cogniterra/Stepik).
    Expects line 1 = k, line 2 = text.
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    k = int(lines[0])
    text = lines[1]

    kmers = kmer_composition(text, k)
    print(f"\nString Composition Result ({len(kmers)} k-mers):")
    for kmer in kmers:
        print(kmer)
    return kmers



if __name__ == "__main__":
    filename = input("Enter input filename (or press Enter for default 'input_1_kmer.txt'): ").strip()
    if not filename:
        filename = "input_1_kmer.txt"

    solve_from_file(filename)