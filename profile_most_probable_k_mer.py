def pr_kmer_profile(kmer, profile):
    """
    Calculates the probability Pr(kmer | Profile) by multiplying 
    the position-specific probabilities for each nucleotide in the k-mer.
    
    profile layout:
    0 or 'A': [prob_0, prob_1, ...]
    1 or 'C': [prob_0, prob_1, ...]
    2 or 'G': [prob_0, prob_1, ...]
    3 or 'T': [prob_0, prob_1, ...]
    """
    nuc_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    prob = 1.0
    
    for i, nuc in enumerate(kmer):
        nuc_row = nuc_map[nuc]
        prob *= profile[nuc_row][i]
        
    return prob


def profile_most_probable_kmer(text, k, profile):
    """
    Finds a Profile-most probable k-mer in Text.
    Ties are broken by returning the first occurring k-mer.
    """
    max_prob = -1.0
    best_kmer = text[:k]  # Default to first k-mer in case all have prob 0

    for i in range(len(text) - k + 1):
        window = text[i:i+k]
        prob = pr_kmer_profile(window, profile)
        
        # Strict '>' ensures we select the FIRST k-mer in case of ties
        if prob > max_prob:
            max_prob = prob
            best_kmer = window

    return best_kmer


def solve_from_file(input_filename="dataset_30305_3.txt", output_filename="output_profile_kmer.txt"):
    """
    Reads dataset file, runs profile_most_probable_kmer, and writes output.
    """
    try:
        with open(input_filename, "r") as f:
            lines = [line.strip() for line in f if line.strip()]

        text = lines[0]
        k = int(lines[1])
        
        # Parse 4xK profile matrix (A, C, G, T)
        profile = []
        for i in range(2, 6):
            profile.append([float(x) for x in lines[i].split()])

        result = profile_most_probable_kmer(text, k, profile)

        with open(output_filename, "w") as f:
            f.write(result + "\n")

        print(f"Result: {result}")
        print(f"Written to: {output_filename}")
        return result
    except FileNotFoundError:
        print(f"Error: File '{input_filename}' not found.")
        return None


def solve_from_string(raw_data_str):
    """
    Parses a raw multi-line dataset string directly without needing a separate file.
    """
    lines = [line.strip() for line in raw_data_str.strip().split("\n") if line.strip()]
    text = lines[0]
    k = int(lines[1])
    
    profile = []
    for i in range(2, 6):
        profile.append([float(x) for x in lines[i].split()])

    result = profile_most_probable_kmer(text, k, profile)
    print(f"Result from string dataset: {result}")
    return result


if __name__ == "__main__":
    # Prompt user for the input filename with a default fallback
    filename = input("Enter input filename (e.g., dataset_30305_3.txt): ").strip()
    if not filename:
        filename = "dataset_30305_3.txt"

    solve_from_file(filename)