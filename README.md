# HW5: Needleman-Wunsch Global Sequence Alignment


## Assignment Overview
The purpose of this assigment is to have you implement the Needleman-Wunsch global pairwise sequence alignment algorithm (dynamic programming).
See this [video](https://www.youtube.com/watch?v=NqYY0PJbD3s) for a walk through of the algorithm implementation. In addition, this is also a helpful resource: https://www.bioinformaticsalgorithms.org/bioinformatics-chapter-5.


## What is Needleman-Wunsch (NW)?
The **Needleman-Wunsch algorithm** is a classic dynamic programming algorithm used in bioinformatics to perform **global pairwise sequence alignment**. It aligns two entire sequences from beginning to end, finding the optimal alignment by considering substitutions, insertions, and deletions (gaps).

### Key Features:
1. **Global Alignment**: Aligns entire sequences from start to finish
2. **Dynamic Programming**: Builds up solutions from smaller subproblems
3. **Gap Penalties**: Uses opening and extension penalties to penalize gap creation and extension
4. **Substitution Matrix**: Uses scoring matrices (like BLOSUM62) to score residue matches/mismatches

### Algorithm Overview:
1. **Initialize three matrices** to track alignment scores:
   - **Alignment matrix**: Scores for matching residues
   - **GapA matrix**: Scores tracking gaps in sequence A
   - **GapB matrix**: Scores tracking gaps in sequence B
2. **Fill matrices forward** using recurrence relations with gap penalties
3. **Backtrace** from the optimal endpoint to reconstruct the actual alignment strings

### Practical Use:
This is essential for comparing DNA/protein sequences across different species to measure evolutionary similarity and relationships.

# Assignment Tasks

## Core Implementation (4 points)
- [ ] **`align/align.py` - `NeedlemanWunsch.align()` method**
  - Initialize the three matrices: alignment matrix, gapA matrix, and gapB matrix
  - Implement the forward-pass dynamic programming algorithm to fill these matrices
  - Use the substitution matrix and gap penalties (open/extend) to calculate scores
  - Return the backtrace result

- [ ] **`align/align.py` - `NeedlemanWunsch._backtrace()` method**
  - Backtrace from the optimal endpoint through the matrices
  - Reconstruct the aligned sequence strings for both seqA and seqB
  - Set the alignment score attribute
  - Return tuple: (alignment_score, seqA_align, seqB_align)

## Main Function (2 points)
- [ ] **`main.py` - Complete `main()` function**
  - Load all 5 species BRD2 sequences from the data folder
  - Align each species to the human BRD2 sequence using:
    - **Substitution matrix**: BLOSUM62
    - **Gap opening penalty**: -10
    - **Gap extension penalty**: -1
  - Print species ordered by alignment score (most similar to least similar)
  - Print the alignment scores for each species-to-human alignment

## Unit Tests (3 points)
- [ ] **`test/test_align.py` - `test_nw_alignment()` function**
  - Load test_seq1.fa and test_seq2.fa
  - Create a NeedlemanWunsch object with BLOSUM62, gap_open=-10, gap_extend=-1
  - Call align() on these sequences
  - Assert that the three matrices are filled correctly
  
- [ ] **`test/test_align.py` - `test_nw_backtrace()` function**
  - Load test_seq3.fa and test_seq4.fa
  - Create a NeedlemanWunsch object with BLOSUM62, gap_open=-10, gap_extend=-1
  - Call align() on these sequences
  - **Expected results:**
    - Alignment score: **17**
    - seqA alignment: `MAVHQLIRRP`
    - seqB alignment: `M---QLIRHP`
  - Assert that the alignments match these expected values

- [ ] **Run all tests with pytest**
  - Execute `pytest test/test_align.py` from the project root
  - Ensure all tests pass

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
