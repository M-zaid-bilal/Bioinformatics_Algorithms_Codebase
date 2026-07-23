# 🧬 Bioinformatics Algorithms — Course 1

**Finding Hidden Messages in DNA**

A hands-on Python implementation of the core algorithms from Week 1 of the *Bioinformatics Specialization*, covering pattern counting, frequent k-mer discovery, reverse complements, clump finding, and genome-scale pattern matching — tested against real bacterial genomes.

---

## 📖 Overview

This repository contains my coursework implementations for **Chapter 1: "Where in the Genome Does DNA Replication Begin?"** — the opening module of the Bioinformatics Specialization. The chapter builds up, algorithm by algorithm, to a real biological question: *where does the origin of replication (`oriC`) hide in a bacterial genome?*

Each script tackles one piece of that puzzle, starting from brute-force string searches and evolving into efficient, dictionary-based algorithms capable of scanning full genomes (including a **4.5 Mb *E. coli*** chromosome) in a fraction of a second. Every file also includes worked examples and standalone test cases so the logic can be verified independently of the course's own grader.

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
│
├── E_coli.txt                             # Full E. coli genome (input data)
├── Vibrio_cholerae.txt                    # Full V. cholerae genome (input data)
├── input_1.txt                            # Sample text/pattern pair
├── input_6.txt                            # Sample text/pattern pair
├── output_positions_vibrio_cholerae.txt   # Saved pattern-matching results
├── output_positions_vibrio_cholerae_2.txt # Saved pattern-matching results
│
├── PatternCount_Chapter_1_Task1.pyproj    # Visual Studio Python project file
├── PatternCount_Chapter_1_Task1.slnx      # Visual Studio solution file
├── .gitattributes
├── .gitignore
└── README.md
```

## 🛠️ Technologies Used

- **Python 3** — standard library only (`os`, `time`, `collections.defaultdict`); no external dependencies
- **Visual Studio** — project developed and managed via `.pyproj` / `.slnx` files (Python workload)
- **Git** — version control

## ⚙️ Setup Instructions

**1. Clone the repository**
```bash
git clone https://github.com/<your-username>/Bioinformatics_Algorithms_Codebase.git
cd Bioinformatics_Algorithms_Codebase
```

**2. Requirements**

Just Python 3 — no `pip install` needed:
```bash
python3 --version
```

**3. Run any script**

Most scripts include a `__main__` block with built-in test cases, so they can simply be run directly:
```bash
python3 Reverse_compliment.py
python3 pattern_matching.py
python3 Frequent_words_efficient.py
python3 Find_Clumps.py
python3 Find_Clumps_Efficient.py
```

> ⚠️ **Note on file paths:** `PatternCount_Chapter_1_Task1.py` and `Find_Clumps_Efficient.py` currently read input files from a hardcoded Windows path (`base_dir = r"C:\Users\ZAID BILAL\..."`). Update `base_dir` to your local repo path (or switch to a relative path like `os.path.dirname(__file__)`) before running these two on another machine.

**4. Import functions into your own scripts**
```python
from Reverse_compliment import reverse_complement
from Frequent_words_efficient import FrequencyTable, BetterFrequentWords
from pattern_matching import pattern_matching, pattern_matching_from_file
from Find_Clumps_Efficient import FindClumpsFast

print(reverse_complement("AAAACCCGGT"))
# → ACCGGGTTTT
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
| Week 2 | Skew Diagrams / Minimum Skew | ⬜ Not started |
| Week 2 | Hamming Distance / Approximate Matching | ⬜ Not started |
| Week 3+ | Motif Finding | ⬜ Not started |

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

## 👤 Author

**Zaid Bilal**
Computer Science, IBA Karachi

## 🙏 Acknowledgments

- **Coursera** and **UC San Diego (UCSD)** for the *Bioinformatics Specialization*
- **Dr. Pavel Pevzner** and **Dr. Phillip Compeau** for the course content and the accompanying textbook *Bioinformatics Algorithms: An Active Learning Approach*
- The **Rosalind**-style genome datasets (*E. coli*, *Vibrio cholerae*) used throughout for real-world testing

## 📄 License

This project is private coursework and is **not licensed for reuse or redistribution**. All rights reserved.