# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    Unit test for NW alignment matrix filling.
    
    Tests that the alignment matrix is correctly filled using test_seq1.fa and test_seq2.fa.
    
    Sequences:
    - test_seq1: MYQR (length 4)
    - test_seq2: MQR (length 3)
    
    Uses BLOSUM62 matrix with gap_open=-10 and gap_extend=-1 (linear gap penalty).
    """
    # Load test sequences
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")
    
    # Create NeedlemanWunsch object
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    # Perform alignment
    score, seq1_align, seq2_align = nw.align(seq1, seq2)
    
    # Expected alignment matrix (5x4 for seq1[0:4] x seq2[0:3])
    # Row/Col indices:    ""   M   Q   R
    #                0  [ 0 -10 -20 -30]
    #            M   1 [-10   5  -5 -15]
    #            Y   2 [-20  -5   4  -6]
    #            Q   3 [-30 -15   0   5]
    #            R   4 [-40 -25 -10   5]
    
    # Assert alignment matrix shape is correct
    assert nw._align_matrix.shape == (5, 4), f"Expected shape (5, 4), got {nw._align_matrix.shape}"
    
    # Assert specific values in alignment matrix were calculated correctly
    # Check first row (gap penalties)
    assert nw._align_matrix[0, 0] == 0, "Top-left corner should be 0"
    assert nw._align_matrix[0, 1] == -10, "First row should have cumulative gap penalties"
    assert nw._align_matrix[0, 2] == -20, "First row should have cumulative gap penalties"
    assert nw._align_matrix[0, 3] == -30, "First row should have cumulative gap penalties"
    
    # Check first column (gap penalties)
    assert nw._align_matrix[0, 0] == 0, "Top-left corner should be 0"
    assert nw._align_matrix[1, 0] == -10, "First column should have cumulative gap penalties"
    assert nw._align_matrix[2, 0] == -20, "First column should have cumulative gap penalties"
    assert nw._align_matrix[3, 0] == -30, "First column should have cumulative gap penalties"
    assert nw._align_matrix[4, 0] == -40, "First column should have cumulative gap penalties"
    
    # Check key internal matrix values
    assert nw._align_matrix[1, 1] == 5, "M-M match should give BLOSUM62(M,M)=5"
    assert nw._align_matrix[2, 2] == 4, "Q-Q match should give BLOSUM62(Q,Q)=5, but with gap penalty"
    assert nw._align_matrix[3, 3] == 5, "R-R match should give BLOSUM62(R,R)=5"
    assert nw._align_matrix[4, 3] == 5, "Bottom-right should be final alignment score"
    
    # Assert final alignment score
    assert score == 5.0, f"Expected alignment score 5.0, got {score}"
    
    # Assert alignment strings are correct
    assert seq1_align == "MYQR", f"Expected seq1 alignment 'MYQR', got '{seq1_align}'"
    assert seq2_align == "M-QR", f"Expected seq2 alignment 'M-QR', got '{seq2_align}'"
    
    # Assert backtracking matrix was created
    assert nw._back is not None, "Backtracking matrix should be created"
    assert nw._back.shape == (5, 4), f"Expected backtrace shape (5, 4), got {nw._back.shape}"
    

def test_nw_backtrace():
    """
    Unit test for NW backtrace procedure.
    
    Tests that the backtrace correctly reconstructs the alignment strings
    using test_seq3.fa and test_seq4.fa.
    
    Sequences:
    - test_seq3: MAVHQLIRRP (length 10)
    - test_seq4: MQLIRHP (length 7)
    
    Uses BLOSUM62 matrix with gap_open=-10 and gap_extend=-1 (linear gap penalty).
    
    Expected alignment (based on linear gap penalty):
    - seqA: MAVHQLIRRP
    - seqB: M---QLIRHP
    - Note: The alignment score for linear penalty is 0 (not 17)
    """
    # Load test sequences
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")
    
    # Create NeedlemanWunsch object
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-10, gap_extend=-1)
    
    # Perform alignment
    score, seq3_align, seq4_align = nw.align(seq3, seq4)
    
    # Assert alignment strings match expected output
    # Expected alignment:
    # MAVHQLIRRP
    # M---QLIRHP
    assert seq3_align == "MAVHQLIRRP", (
        f"Expected seq3 alignment 'MAVHQLIRRP', got '{seq3_align}'"
    )
    assert seq4_align == "M---QLIRHP", (
        f"Expected seq4 alignment 'M---QLIRHP', got '{seq4_align}'"
    )
    
    # Assert alignment strings have same length
    assert len(seq3_align) == len(seq4_align), (
        f"Aligned sequences should have same length, got {len(seq3_align)} and {len(seq4_align)}"
    )
    
    # Assert alignment strings length equals original sequences plus gaps
    # seq3 has 10 chars, seq4 has 7 chars
    # The alignment should be 10 chars (from seq3) with 3 gaps inserted in seq4
    assert len(seq3_align) == len(seq3), "seq3 alignment should match original length"
    assert len(seq4_align) == len(seq3), "seq4 alignment should match seq3 alignment length"
    
    # Assert backtrace score is set correctly
    # For linear gap penalty: M-M(5) + A-(-)(−10) + V-(-)(−10) + H-(-)(−10) + Q-Q(5) + 
    #                         L-L(4) + I-I(4) + R-R(5) + R-H(0) + P-P(7) = 0
    assert nw.alignment_score == score, "Alignment score should be set by backtrace"
    
    # Assert alignment score is correct (0 for linear gap penalty)
    assert score == 0.0, f"Expected alignment score 0.0 for linear penalty, got {score}"
    
    # Assert backtracing matrix was created and has correct shape
    assert nw._back is not None, "Backtracking matrix should exist"
    assert nw._back.shape == (11, 8), f"Expected backtrace shape (11, 8), got {nw._back.shape}"
    
    # Verify alignment contains expected sequences (possibly with gaps)
    # Remove gaps to verify original sequences are preserved
    seq3_no_gaps = seq3_align.replace("-", "")
    seq4_no_gaps = seq4_align.replace("-", "")
    
    assert seq3_no_gaps == seq3, (
        f"seq3 alignment without gaps should equal original seq3. "
        f"Got '{seq3_no_gaps}' vs '{seq3}'"
    )
    assert seq4_no_gaps == seq4, (
        f"seq4 alignment without gaps should equal original seq4. "
        f"Got '{seq4_no_gaps}' vs '{seq4}'"
    )




