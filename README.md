# HW5: Needleman-Wunsch Global Sequence Alignment

For HW5, please ignore the references to an affine gap penalty in the assignment. Implement the NW algorithm with a linear gap penalty. If you’ve already implemented the affine penalty version we’ll give a little EC for it. We will update the test case score. Sorry for this oversight.

## Assignment Overview
The purpose of this assigment is to have you implement the Needleman-Wunsch global pairwise sequence alignment algorithm (dynamic programming).
See this [video](https://www.youtube.com/watch?v=NqYY0PJbD3s) for a walk through of the algorithm implementation. In addition, this is also a helpful resource: https://www.bioinformaticsalgorithms.org/bioinformatics-chapter-5.


## What is Needleman-Wunsch (NW)?
The **Needleman-Wunsch algorithm** is a classic dynamic programming algorithm used in bioinformatics to perform **global pairwise sequence alignment**. It aligns two entire sequences from beginning to end, finding the optimal alignment by considering substitutions, insertions, and deletions (gaps).

### Key Features:
1. **Global Alignment**: Aligns entire sequences from start to finish
2. **Dynamic Programming**: Builds up solutions from smaller subproblems
3. **Linear Gap Penalty**: Uses a single penalty for any gap (no distinction between opening and extending)
4. **Substitution Matrix**: Uses scoring matrices (like BLOSUM62) to score residue matches/mismatches

### Algorithm Overview:
1. **Initialize alignment matrix** to track alignment scores:
   - **Alignment matrix**: Scores for matching residues and gaps combined
2. **Fill matrix forward** using recurrence relations with linear gap penalty (gap_cost = gap_penalty × gap_length)
3. **Backtrace** from the optimal endpoint to reconstruct the actual alignment strings

### Important Note on Gap Penalty:
 **Use LINEAR gap penalty** (not affine). This means:
- Any gap of length *n* costs: `n × gap_penalty`
- Do NOT distinguish between gap opening and gap extension
- The `gap_extend` parameter in the class is ignored for linear penalty implementation

### Practical Use:
This is essential for comparing DNA/protein sequences across different species to measure evolutionary similarity and relationships.

# Assignment Tasks

## 1. Core Implementation (4 points)
- [x] **`align/align.py` - `NeedlemanWunsch.align()` method**
  - Initialize the alignment matrix (single matrix, NOT three separate matrices)
  - Implement the forward-pass dynamic programming algorithm with **LINEAR gap penalty**
  - Use the substitution matrix and gap penalty to calculate scores
  - Gap cost = gap_penalty × gap_length (e.g., if gap_open=-10, a 2-char gap costs -20)
  - Return the backtrace result

- [x] **`align/align.py` - `NeedlemanWunsch._backtrace()` method**
  - Backtrace from the optimal endpoint through the matrix
  - Reconstruct the aligned sequence strings for both seqA and seqB
  - Set the alignment score attribute
  - Return tuple: (alignment_score, seqA_align, seqB_align)



```
cd /Users/jugomez/Library/CloudStorage/Box-Box/Biological_and_Medical_Informatics/Biocomputing_Algorithms_\(BMI_203\)/Assignments/HW5-NW && python -c "
from align import NeedlemanWunsch, read_fasta

# Load test sequences
seq3, _ = read_fasta('./data/test_seq3.fa')
seq4, _ = read_fasta('./data/test_seq4.fa')

# Create NW object with BLOSUM62 and linear gap penalty
nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-4, gap_extend=-1)

# Align sequences
score, seqA_align, seqB_align = nw.align(seq3, seq4)

print(f'Alignment Score: {score}')
print(f'seqA: {seqA_align}')
print(f'seqB: {seqB_align}')
print()
print('Expected:')
print(f'Score: 1')
print(f'seqA: MAVHQLIRRP')
print(f'seqB: M---QLIRHP')
"
```8

```
Alignment Score: 18
seqA: MAVHQLIRRP
seqB: M---QLIRHP

Expected:
Score: 18
seqA: MAVHQLIRRP
seqB: M

```

```
M-M:     +5
A-gap:  -10
V-gap:  -10
H-gap:  -10
Q-Q:     +5
L-L:     +4
I-I:     +4
R-R:     +5
R-H:      0
P-P:     +7
───────────
Total:    0
```

## 2. Main Function (2 points)
- [x] **`main.py` - Complete `main()` function**
  - Load all 5 species BRD2 sequences from the data folder
  - Align each species to the human BRD2 sequence using:
    - **Substitution matrix**: BLOSUM62
    - **Gap penalty**: -10 (linear - use only the gap_open parameter)
  - Print species ordered by alignment score (most similar to least similar)
  - Print the alignment scores for each species-to-human alignment


```
cd /Users/jugomez/Library/CloudStorage/Box-Box/Biological_and_Medical_Informatics/Biocomputing_Algorithms_\(BMI_203\)/Assignments/HW5-NW && python main.py
```

* Dolphin is most evolutionarily similar to humans (mammals)

* Mouse is also a mammal, so high similarity makes sense

* Chicken (bird) has lower similarity than mammals

* Shoebill (bird) has the lowest similarity among tested species


## 3. Unit Tests (3 points)
- [x] **`test/test_align.py` - `test_nw_alignment()` function**
  - Load test_seq1.fa and test_seq2.fa
  - Create a NeedlemanWunsch object with BLOSUM62, gap_open=-10
  - Call align() on these sequences
  - Assert that the alignment matrix is filled correctly
  
- [x] **`test/test_align.py` - `test_nw_backtrace()` function**
  - Load test_seq3.fa and test_seq4.fa
  - Create a NeedlemanWunsch object with BLOSUM62, gap_open=-10
  - Call align() on these sequences
  - **Expected results:**
    - Alignment score: **17**
    - seqA alignment: `MAVHQLIRRP`
    - seqB alignment: `M---QLIRHP`
  - Assert that the alignments match these expected values

- [x] **Run all tests with pytest**
  - Execute `pytest test/test_align.py` from the project root
  - Ensure all tests pass


```
cd /Users/jugomez/Library/CloudStorage/Box-Box/Biological_and_Medical_Informatics/Biocomputing_Algorithms_\(BMI_203\)/Assignments/HW5-NW && python -c "
from align import NeedlemanWunsch, read_fasta
import numpy as np

# Test seq1 and seq2
seq1, _ = read_fasta('./data/test_seq1.fa')
seq2, _ = read_fasta('./data/test_seq2.fa')

print(f'seq1: {seq1}')
print(f'seq2: {seq2}')
print()

nw = NeedlemanWunsch('./substitution_matrices/BLOSUM62.mat', gap_open=-10, gap_extend=-1)
score, seq1_align, seq2_align = nw.align(seq1, seq2)

print(f'Alignment score: {score}')
print(f'seq1 alignment: {seq1_align}')
print(f'seq2 alignment: {seq2_align}')
print()
print('Alignment matrix shape:', nw._align_matrix.shape)
print('Alignment matrix:')
print(nw._align_matrix)
"

```

```
cd /Users/jugomez/Library/CloudStorage/Box-Box/Biological_and_Medical_Informatics/Biocomputing_Algorithms_\(BMI_203\)/Assignments/HW5-NW && "/Users/jugomez/Library/CloudStorage/Box-Box/Biological_and_Medical_Informatics/Biocomputing_Algorithms_(BMI_203)/Assignments/HW5-NW/.venv/bin/python" -m pytest test/test_align.py -v
```



## Package Setup (1 point)
- [ ] **Create `pyproject.toml` configuration file**
  - Use flit as the build backend
  - Make the package installable with `pip install .`
  - Include necessary metadata and dependencies

## Code Style (1 point)
- [ ] **Ensure code quality**
  - Add clear comments explaining algorithm logic
  - Include docstrings for all methods
  - Follow PEP 8 style guidelines
  - Make code readable and maintainable

## Extra Credit (0.5 points)
- [ ] **GitHub Actions/Workflow setup** (optional)
  - Create CI/CD pipeline to automatically run tests on push

# Getting Started
To get started you will need to fork this repository onto your own Github account. Work on the codebase from your own repo and commit changes. 

The following packages will be needed:
* numpy
* pytest

# Completing the assignment
Make sure to push all your code to Github, ensure that your unit tests are correct, and submit a link to your Github through the Google classroom assignment.

# Grading
## Code (6 points)
* Pairwise global alignment works properly (6)
    * Correct implementation of Needleman-Wunsch algorithm (4)
    * Produces correct order of species in main.py (1) 
    * Produces correct NW alignment scores in main.py (1)

## Unit tests (3 points)
* `test_nw_alignment` function properly checks that matrices are filled in correctly for alignment of test_seq1.fa and test_seq2.fa (1)
* `test_nw_backtrace` function properly checks that backtrace works correctly (1)
* Ensure functionality with pytest (1)
## Style (1 points)
* Readable code with clear comments and method descriptions (1)
## Extra credit (0.5)
* Github actions/workflow (0.5)

























