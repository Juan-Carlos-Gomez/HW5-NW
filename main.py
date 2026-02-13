# Import NeedlemanWunsch class and read_fasta function
from align import read_fasta, NeedlemanWunsch

def main():
    """
    This function aligns all species BRD2 sequences to the human BRD2 sequence
    and prints them in order of similarity (most to least similar).
    
    It also prints the alignment scores for each species-to-human alignment.
    
    Uses:
    - Substitution matrix: BLOSUM62
    - Gap penalty: -10 (linear gap penalty)
    """
    # Load all species BRD2 sequences
    hs_seq, hs_header = read_fasta("./data/Homo_sapiens_BRD2.fa")
    gg_seq, gg_header = read_fasta("./data/Gallus_gallus_BRD2.fa")
    mm_seq, mm_header = read_fasta("./data/Mus_musculus_BRD2.fa")
    br_seq, br_header = read_fasta("./data/Balaeniceps_rex_BRD2.fa")
    tt_seq, tt_header = read_fasta("./data/Tursiops_truncatus_BRD2.fa")

    # Create a list to store species info: (species_name, sequence, header)
    species_list = [
        ("Gallus gallus (Chicken)", gg_seq, gg_header),
        ("Mus musculus (Mouse)", mm_seq, mm_header),
        ("Balaeniceps rex (Shoebill)", br_seq, br_header),
        ("Tursiops truncatus (Dolphin)", tt_seq, tt_header),
    ]
    
    # Create NeedlemanWunsch object with BLOSUM62 and linear gap penalty
    nw = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat", gap_open=-4, gap_extend=-4)
    
    # Align each species to human and store results
    alignment_results = []
    
    for species_name, species_seq, species_header in species_list:
        # Align species to human
        score, hs_align, species_align = nw.align(hs_seq, species_seq)
        
        # Store results as tuple: (score, species_name, human_alignment, species_alignment)
        alignment_results.append((score, species_name, hs_align, species_align))
    
    # Sort by alignment score in descending order (most similar first)
    alignment_results.sort(key=lambda x: x[0], reverse=True)
    
    # Print header
    print("=" * 80)
    print("Needleman-Wunsch Global Alignment: Species BRD2 to Human BRD2")
    print("=" * 80)
    print()
    
    # Print species in order of similarity (most to least similar)
    print("Species ranked by similarity to Human BRD2:")
    print("-" * 80)
    for rank, (score, species_name, _, _) in enumerate(alignment_results, 1):
        print(f"{rank}. {species_name:40} | Score: {score:8.1f}")
    
    print()
    print("=" * 80)
    print("Detailed Alignments")
    print("=" * 80)
    
    # Print detailed alignment for each species
    for score, species_name, hs_align, species_align in alignment_results:
        print()
        print(f"\n{species_name}")
        print("-" * 80)
        print(f"Alignment Score: {score}")
        print()
        print("Human BRD2:")
        print(hs_align)
        print()
        print(f"{species_name}:")
        print(species_align)
        print()
    

if __name__ == "__main__":
    main()
