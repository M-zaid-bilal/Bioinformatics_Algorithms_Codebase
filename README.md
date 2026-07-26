# 🧬 Bioinformatics Algorithms — Course 1

**Finding Hidden Messages in DNA**

A hands-on Python implementation of the core algorithms from Weeks 1–2 of the *Bioinformatics Specialization*, covering pattern counting, frequent k-mer discovery, reverse complements, clump finding, GC skew analysis, and approximate (mismatch-tolerant) pattern matching — culminating in a full DnaA-box finding pipeline tested against real bacterial genomes.

---

## 📖 Overview

This repository contains my coursework implementations for **Chapter 1: "Where in the Genome Does DNA Replication Begin?"** — the opening module of the Bioinformatics Specialization. The chapter builds up, algorithm by algorithm, to a real biological question: *where does the origin of replication (`oriC`) hide in a bacterial genome?*

Each script tackles one piece of that puzzle, starting from brute-force string searches and evolving into efficient, dictionary-based algorithms capable of scanning full genomes (including a **4.5 Mb *E. coli*** chromosome) in a fraction of a second. The project now goes past exact-match searching into **approximate pattern matching with mismatches**, **GC skew analysis**, and finally combines all of it into `dna_finder.py` — an end-to-end pipeline that locates the likely `oriC` region and candidate DnaA boxes directly from a raw genome file, complete with a generated skew-diagram plot. Every file also includes worked examples and standalone test cases so the logic can be verified independently of the course's own grader.

## 🎓 Course Information

| | |
|---|---|
| **Specialization** | [Bioinformatics Specialization](https://www.coursera.org/specializations/bioinformatics) |
| **Course** | Course 1 — Finding Hidden Messages in DNA (Bioinformatics I) |
| **University** | University of California, San Diego (UCSD) |
| **Platform** | Coursera |
| **Instructors** | Pavel Pevzner, Phillip Compeau |

## 🧩 Topics Covered

| Algorithm | File | Description |
|---|---|---|
| **Pattern Count** | `PatternCount_Chapter_1_Task1.py` | Brute-force count of how many times a pattern occurs (with overlaps) in a text |
| **Frequent Words (naive)** | `Frequent_words.py` | Finds the most frequent *k*-mers in a string by counting every substring against every other |
| **Frequency Table / Better Frequent Words** | `Frequent_words_efficient.py` | Builds a hash-map frequency table in one pass, then extracts the max-count *k*-mers — the efficient rewrite of the naive version above |
| **Reverse Complement** | `Reverse_compliment.py` | Reverses a DNA strand and swaps each base for its Watson–Crick complement (A↔T, C↔G) |
| **Pattern Matching** | `pattern_matching.py` | Returns every starting position of a pattern in a genome, including a FASTA-aware file reader for genome-scale search |
| **Clump Finding (naive)** | `Find_Clumps.py` | Slides an L-length window across the genome and flags *k*-mers appearing ≥ *t* times inside it — an early flag for `oriC` |
| **Clump Finding (efficient)** | `Find_Clumps_Efficient.py` | Pre-indexes every *k*-mer's positions once, then checks clump membership in a single pass — fast enough for the full *E. coli* genome |
| **GC Skew / Minimum Skew** | `Skew_Count.py` | Tracks the running G−C count across a genome and returns the position(s) of minimum skew — the classic first estimate of `oriC`'s location |
| **Hamming Distance / Approximate Pattern Matching** | `Hamming_Distance.py` | Counts mismatched positions between equal-length strings, then uses it to find every position in a genome where a pattern occurs with ≤ *d* mismatches |
| **d-Neighborhood Generation** | `Immediate_neighbor.py` | Recursively generates every *k*-mer within Hamming distance *d* of a given pattern (the neighborhood used by the mismatch-tolerant searches below) |
| **Frequent Words with Mismatches** | `Freequent_words_with_mismatches.py` | Extends Better Frequent Words to tolerate up to *d* mismatches per occurrence, using the neighborhood generator + approximate counting |
| **Frequent Words with Mismatches + Reverse Complement** | `Frequent_Words_with_mismatches_and_reverse_complement.py` | Same as above, but also counts each k-mer's reverse-complement neighborhood — the exact motif definition of a biological DnaA box |
| **DnaA Box / oriC Finder Pipeline** | `dna_finder.py` | End-to-end script: reads a raw genome file → computes minimum skew to estimate `oriC` → searches a window around it for mismatch- and RC-tolerant DnaA box candidates → maps their genome positions → plots and saves a GC-skew diagram (`gc_skew_plot.png`) → writes a full JSON report (`dnaa_results.json`) |

## 🌳 File Structure

```
Bioinformatics_Algorithms_Codebase-master/
├── PatternCount_Chapter_1_Task1.py        # Core pattern_count() + file-based test runner
├── Frequent_words.py                      # Naive FrequentWords algorithm
├── Frequent_words_efficient.py            # FrequencyTable + BetterFrequentWords
├── Reverse_compliment.py                  # reverse_complement()
├── pattern_matching.py                    # pattern_matching() + FASTA file search
├── Find_Clumps.py                         # Naive (L, t)-clump finder
├── Find_Clumps_Efficient.py               # Optimized clump finder (genome-scale)
├── Skew_Count.py                          # calculate_skew() + minimum_skew()
├── Hamming_Distance.py                    # hamming_distance() + approximate_pattern_matching()
├── Immediate_neighbor.py                  # neighbors() — d-neighborhood generator
├── Freequent_words_with_mismatches.py     # approximate_pattern_count() + mismatch-tolerant frequent words
├── Frequent_Words_with_mismatches_and_reverse_complement.py  # + reverse-complement-aware variant
├── dna_finder.py                          # Full oriC/DnaA-box pipeline (skew + motif search + plotting)
│
├── E_coli.txt                             # Full E. coli genome (input data)
├── Vibrio_cholerae.txt                    # Full V. cholerae genome (input data)
├── Salmonella_enterica.txt                # Full S. enterica genome (input data, used by dna_finder.py)
├── input_1.txt, input_6.txt               # Sample text/pattern pairs
├── input_1_app.txt, input_1_approximate_pattern_count.txt,
│   input_1_fwm.txt, input_1_fwm_rc.txt    # Sample inputs for approximate-matching/mismatch scripts
├── dataset_30278_*.txt                    # Course grader datasets (Hamming distance, neighbors, etc.)
├── output_positions_vibrio_cholerae.txt   # Saved pattern-matching results
├── output_positions_vibrio_cholerae_2.txt # Saved pattern-matching results
├── dnaa_results.json                      # Saved output of dna_finder.py (DnaA box candidates + positions)
├── gc_skew_plot.png                       # Saved GC-skew diagram generated by dna_finder.py
│
├── PatternCount_Chapter_1_Task1.pyproj    # Visual Studio Python project file
├── PatternCount_Chapter_1_Task1.slnx      # Visual Studio solution file
├── .gitattributes
├── .gitignore
└── README.md
```

## 🛠️ Technologies Used

- **Python 3** — mostly standard library (`os`, `time`, `json`, `collections.defaultdict`)
- **Matplotlib** — used by `dna_finder.py` to render and save the GC-skew diagram (`pip install matplotlib`)
- **Visual Studio** — project developed and managed via `.pyproj` / `.slnx` files (Python workload)
- **Git** — version control

## ⚙️ Setup Instructions

**1. Clone the repository**
```bash
git clone https://github.com/<your-username>/Bioinformatics_Algorithms_Codebase.git
cd Bioinformatics_Algorithms_Codebase
```

**2. Requirements**

Python 3 for everything; `matplotlib` is only needed if you want to run `dna_finder.py` (it generates the GC-skew plot):
```bash
python3 --version
pip install matplotlib   # only required for dna_finder.py
```

**3. Run any script**

Most scripts include a `__main__` block with built-in test cases, so they can simply be run directly:
```bash
python3 Reverse_compliment.py
python3 pattern_matching.py
python3 Frequent_words_efficient.py
python3 Find_Clumps.py
python3 Find_Clumps_Efficient.py
python3 Skew_Count.py
python3 Hamming_Distance.py
python3 Immediate_neighbor.py
python3 Freequent_words_with_mismatches.py
python3 Frequent_Words_with_mismatches_and_reverse_complement.py

# Full oriC / DnaA-box finder pipeline (needs matplotlib):
pip install matplotlib
python3 dna_finder.py
```

> ⚠️ **Note on file paths:** `PatternCount_Chapter_1_Task1.py` and `Find_Clumps_Efficient.py` currently read input files from a hardcoded Windows path (`base_dir = r"C:\Users\ZAID BILAL\..."`). Update `base_dir` to your local repo path (or switch to a relative path like `os.path.dirname(__file__)`) before running these two on another machine.

**4. Import functions into your own scripts**
```python
from Reverse_compliment import reverse_complement
from Frequent_words_efficient import FrequencyTable, BetterFrequentWords
from pattern_matching import pattern_matching, pattern_matching_from_file
from Find_Clumps_Efficient import FindClumpsFast
from Skew_Count import calculate_skew, minimum_skew
from Hamming_Distance import hamming_distance, approximate_pattern_matching
from Immediate_neighbor import neighbors
from Frequent_Words_with_mismatches_and_reverse_complement import frequent_words_with_mismatches_and_rc
from dna_finder import find_dna_a_box

print(reverse_complement("AAAACCCGGT"))
# → ACCGGGTTTT

# Run the full pipeline on any genome file:
results = find_dna_a_box("Salmonella_enterica.txt", k=9, d=1, window_size=500)
```

## ✅ Progress Tracker

| Week | Topic | Status |
|---|---|---|
| Week 1 | Pattern Count | ✅ Complete |
| Week 1 | Frequent Words (naive) | ✅ Complete |
| Week 1 | Frequency Table / Better Frequent Words | ✅ Complete |
| Week 1 | Reverse Complement | ✅ Complete |
| Week 1 | Pattern Matching (incl. genome file search) | ✅ Complete |
| Week 1 | Clump Finding (naive) | ✅ Complete |
| Week 1 | Clump Finding (efficient, genome-scale) | ✅ Complete |
| Week 2 | Skew Diagrams / Minimum Skew | ✅ Complete |
| Week 2 | Hamming Distance / Approximate Pattern Matching | ✅ Complete |
| Week 2 | d-Neighborhood Generation | ✅ Complete |
| Week 2 | Frequent Words with Mismatches (+ Reverse Complement) | ✅ Complete |
| Week 2 | DnaA Box / oriC Finder Pipeline (skew + motif search + plotting) | ✅ Complete — tested end-to-end on *Salmonella enterica* |
| Week 3 | Motif Enumeration (Brute-Force Motif Search) | ⬜ Not started |
| Week 3 | Median String / Profile-most Probable k-mer | ⬜ Not started |
| Week 3 | Greedy Motif Search (with Laplace/pseudocount smoothing) | ⬜ Not started |
| Week 4 | Randomized Motif Search | ⬜ Not started |
| Week 4 | Gibbs Sampling | ⬜ Not started |
| Week 5 | Bioinformatics Application Challenge (real motif dataset, *M. tuberculosis*) | ⬜ Not started |

## 💻 Sample Outputs

**Pattern Count**
```python
pattern_count("ATGCATTCATCG", "TC")
# → 2
```

**Reverse Complement**
```python
reverse_complement("AAAACCCGGT")
# → "ACCGGGTTTT"
```

**Better Frequent Words**
```python
BetterFrequentWords("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4)
# → "GCAT CATG"
```

**Pattern Matching**
```python
pattern_matching("GATATATGCATATACTT", "ATAT")
# → [1, 3, 9]
```

**Genome-Scale Pattern Matching** (searching *Vibrio cholerae* for the `oriC`-associated 9-mer `CTTGATCAT`)
```python
pattern_matching_from_file("Vibrio_cholerae.txt", "CTTGATCAT")
# → [116556, 149355, 151913, 152013, 152394, 186189, 194276,
#    200076, 224527, 307692, 479770, 610980, 653338, 679985,
#    768828, 878903, 985368]
```
Notice how several matches cluster tightly around position ~152,000 — exactly the kind of clump that points to a replication origin.

**Efficient Clump Finding on the full *E. coli* genome** (*k*=9, L=500, t=3)
```
Number of distinct k-mers forming (L, t)-clumps in E. coli genome: 1904
Time taken to find clumps: <1 second
```

**Minimum Skew**
```python
minimum_skew("GAGCCACCGCGATA")
# → [9, 10, 11, 12]
```

**Approximate Pattern Matching (Hamming distance ≤ d)**
```python
approximate_pattern_matching("CGCCCGAATCCAGAACGCATTCCCATGTACA", "ATTCTGGA", 3)
# → positions where "ATTCTGGA" occurs with at most 3 mismatches
```

**d-Neighborhood Generation**
```python
neighbors("GGGGGGCGT", 2)
# → set of every 9-mer within Hamming distance 2 of "GGGGGGCGT"
```

**Frequent Words with Mismatches + Reverse Complement**
```python
frequent_words_with_mismatches_and_rc("ACGTTGCATGTCGCATGATGCATGAGAGCT", 4, 1)
# → most frequent 4-mers counting both a k-mer's mismatch-neighborhood
#   and its reverse complement's — the actual definition of a DnaA box motif
```

**DnaA Box / oriC Finder Pipeline** (*Salmonella enterica*, *k*=9, *d*=1, window=500)
```
Genome Length: 4,809,037 bp
Estimated oriC Center: 3,764,856
Window Analyzed: 3,764,606 to 3,765,106
Candidate DnaA Boxes: ['CCGGAAGCT', 'AGCTTCCGG']
  > Motif 'CCGGAAGCT': [3764688, 3764693, 3764747, 3764998, 3765003]
  > Motif 'AGCTTCCGG': [3764688, 3764693, 3764747, 3764998, 3765003]
```
Note that `AGCTTCCGG` is the reverse complement of `CCGGAAGCT` — both mapping to the exact same five genome positions confirms the pipeline is correctly picking up a real palindromic-ish DnaA box cluster right at the minimum-skew point. Results are saved to `dnaa_results.json`, and the full genome-wide skew curve (with the oriC window highlighted) is saved to `gc_skew_plot.png`.

## 🔭 Next Steps

Having finished Chapter 1 in full (Weeks 1–2: pattern warm-up algorithms → skew analysis → mismatch-tolerant motif search → a working `oriC`/DnaA-box pipeline), the natural continuation is **Chapter 2: "Which DNA Patterns Play the Role of Molecular Clocks?"**, which shifts from replication origins to regulatory motif discovery:

- **Week 3 — Hunting for Regulatory Motifs:** Brute-force Motif Enumeration, the Median String problem, and Profile-most Probable k-mer / Greedy Motif Search (with Laplace pseudocount smoothing) for finding shared motifs across a set of DNA strings rather than within one genome.
- **Week 4 — Randomized Algorithms:** Randomized Motif Search and Gibbs Sampling, which "roll dice" to converge on motifs far faster than the deterministic Week 3 methods scale.
- **Week 5 — Bioinformatics Application Challenge:** applying motif-finding software to a real dataset — genes thought to help *Mycobacterium tuberculosis* enter dormancy in a host.

Planned repo additions to match: `Motif_Enumeration.py`, `Median_String.py`, `Greedy_Motif_Search.py`, `Randomized_Motif_Search.py`, and `Gibbs_Sampler.py`, each following the same pattern used throughout this repo — a pure algorithm function, a file-based test runner, and inline worked examples.

## 👤 Author

**Zaid Bilal**
Computer Science, IBA Karachi

## 🙏 Acknowledgments

- **Coursera** and **UC San Diego (UCSD)** for the *Bioinformatics Specialization*
- **Dr. Pavel Pevzner** and **Dr. Phillip Compeau** for the course content and the accompanying textbook *Bioinformatics Algorithms: An Active Learning Approach*
- The **Rosalind**-style genome datasets (*E. coli*, *Vibrio cholerae*, *Salmonella enterica*) used throughout for real-world testing

## 📄 License

This project is private coursework and is **not licensed for reuse or redistribution**. All rights reserved.