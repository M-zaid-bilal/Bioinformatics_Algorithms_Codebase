from collections import defaultdict
import Frequent_words_efficient as fwe
import os 
import time
def FindClumpsFast(Text, k, L, t):
    patterns = set()
    positions = defaultdict(list)
    
    # Store starting positions of all k-mers
    for i in range(len(Text) - k + 1):
        kmer = Text[i : i + k]
        positions[kmer].append(i)
        
    # Check only candidates that appear at least t times in the whole genome
    for kmer, pos_list in positions.items():
        if len(pos_list) >= t:
            for i in range(len(pos_list) - t + 1):
                # Check if t occurrences fit inside window of size L
                if (pos_list[i + t - 1] + k - pos_list[i]) <= L:
                    patterns.add(kmer)
                    break
                    
    return patterns

#Now let's test it on E.coli genome data
if __name__ == "__main__":
    # Example usage with E.coli genome data
    base_dir = r"C:\Users\ZAID BILAL\source\repos\PatternCount_Chapter_1_Task1"
    input_file = "E_coli.txt"  # Replace with your actual file name
    filepath = os.path.join(base_dir, input_file)
    #Error handling if Ecoli is not present

    with open(filepath, "r") as f:
        Text = f.read().strip()
    
    k = 9  # Length of k-mer
    L = 500  # Length of the window
    t = 3   # Minimum number of occurrences
    #FIND the time taken to find clumps in E.coli genome
 
    start_time = time.time()
    clumps = FindClumpsFast(Text, k, L, t)
    end_time = time.time()
    no_of_clumps= len(clumps)
    #instead of printing whole set, let's print the count of distinct k-mers forming (L, t)-clumps
    print(f"Number of distinct k-mers forming (L, t)-clumps in E.coli genome:", no_of_clumps)
    print(f"Time taken to find clumps: {end_time - start_time:.4f} seconds")