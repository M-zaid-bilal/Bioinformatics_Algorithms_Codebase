from Freequent_words_with_mismatches import neighbors
from Hamming_Distance import hamming_distance

def motif_enumeration(dna, k, d):
    """Finds all (k, d)-motifs present in all DNA sequences."""
    patterns = set()
    first_string = dna[0]
    
    # 1. Generate candidate k-mers and their neighbors from the FIRST DNA sequence
    candidate_kmers = set()
    for i in range(len(first_string) - k + 1):
        kmer = first_string[i:i+k]
        candidate_kmers.update(neighbors(kmer, d))
    
    # 2. Test each candidate against ALL DNA sequences
    for candidate in candidate_kmers:
        is_motif_in_all = True
        
        for string in dna:
            found = False
            for i in range(len(string) - k + 1):
                sub_kmer = string[i:i+k]
                if hamming_distance(candidate, sub_kmer) <= d:
                    found = True
                    break
            
            if not found:
                is_motif_in_all = False
                break
                
        if is_motif_in_all:
            patterns.add(candidate)
            
    return sorted(list(patterns))

def get_input_interactively():
    """Prompts the user to enter k, d, and DNA strings manually in terminal."""
    print("\n--- Manual Input Mode ---")
    k = int(input("Enter k (length of motif): ").strip())
    d = int(input("Enter d (max mismatches): ").strip())
    
    print("Enter DNA sequences separated by space or press Enter after each sequence.")
    print("Type 'done' on a new line when you are finished entering sequences:")
    
    dna = []
    while True:
        line = input().strip()
        if line.lower() == 'done':
            break
        if line:
            # Handle if multiple space-separated sequences were pasted on one line
            dna.extend(line.split())
            
    return k, d, dna

def get_input_from_file(filename="input_6_motif_enum.txt"):
    """Reads k, d, and DNA strings from a file."""
    with open(filename, "r") as f:
        tokens = f.read().strip().split()
    
    if not tokens:
        raise ValueError(f"{filename} is empty!")
        
    k = int(tokens[0])
    d = int(tokens[1])
    dna = tokens[2:]
    return k, d, dna

def main():
    print("Choose input method:")
    print("1. Read from dataset.txt")
    print("2. Enter manually in terminal")
    
    choice = input("Choice (1 or 2): ").strip()
    
    if choice == "2":
        k, d, dna = get_input_interactively()
    else:
        try:
            k, d, dna = get_input_from_file("dataset.txt")
        except FileNotFoundError:
            print("dataset.txt not found! Switching to manual terminal input...")
            k, d, dna = get_input_interactively()

    if not dna:
        print("Error: No DNA sequences provided!")
        return

    print(f"\nProcessing: k={k}, d={d}, Sequences={len(dna)}")
    
    result_motifs = motif_enumeration(dna, k, d)
    output_string = " ".join(result_motifs)
    
    print("\nResult:")
    print(output_string)
    
    # Save output to file
    with open("output.txt", "w") as f:
        f.write(output_string)
    print("\n(Result also saved to output.txt)")

if __name__ == "__main__":
    main()