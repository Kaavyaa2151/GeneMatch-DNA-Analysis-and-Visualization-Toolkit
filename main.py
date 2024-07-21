# Loading the Sequence Class from the "Sequence.py"
from Sequence import Sequence
from print_string_sequence import print_string_sequence
from sw_algo import sw_algo
from bm_approx_mat import bm_approx_mat
from detmut import detmut
from create_phylogenic_tree import create_phylogenic_tree

# Defining the Hash Table for Smith-Waterman alignment
# "N" is the value of the neuclotide that has been incorrectly identified
score_table_sw = {
    "A": {"A": 2, "C": -1, "G": -1, "T": -1, "N": -1},
    "C": {"A": -1, "C": 2, "G": -1, "T": -1, "N": -1},
    "G": {"A": -1, "C": -1, "G": 2, "T": -1, "N": -1},
    "T": {"A": -1, "C": -1, "G": -1, "T": 2, "N": -1},
    "N": {"A": -1, "C": -1, "G": -1, "T": -1, "N": -1},
    "gap": {"gap": -1}
}

# Loading two sequences from data files
seq1 = Sequence("data/human_seq1.fasta")
seq2 = Sequence("data/human_seq2.fasta")

# Printing the size of the sequences
print("Size of Sequence 1 is:", seq1.getSize())
print("Size of Sequence 2 is:", seq2.getSize())

# Smith-Waterman alignment
alseqs = sw_algo(seq1, seq2, score_table_sw)

# Printing Results
print("\n\n")
print("Smith-Waterman Alignment:")
print("===============================================================")
print("Sequence 1:")
print_string_sequence(alseqs[0])
print("\nSequence 2:")
print_string_sequence(alseqs[0])
print("===============================================================")

# Loading a gene pattern to be identified
pattern = Sequence("data/pattern.fasta")

# Boyer-Moore approximate matching
approx_match_positions = bm_approx_mat(pattern, seq1, maxmis=2)

# Printing Results
print("\n\n")
print("Boyer-Moore Approximate Matches:")
print("===============================================================")
if approx_match_positions:
    for pos in approx_match_positions:
        print(f"Approximate match starting at position {pos}")
else:
    print("No Boyer-Moore approximate matches detected.")
print("===============================================================\n")

# Mutation detection using Hamming distance
# Printing Results
print("\n\n")
print("Hamming Distance for Mutation Detection:")
print("===============================================================")
detmut(seq1, seq2)
print("===============================================================\n")

# Create Phylogenetic Tree
species = {
    "human": Sequence("data/human_seq1.fasta"),
    "frog": Sequence("data/blue_poison_arrow_frog.fasta"),
    "elephant": Sequence("data/african_elephant.fasta"),
    "whale": Sequence("data/blue_whale.fasta"),
    "dog": Sequence("data/dog.fasta"),
    "pigeon": Sequence("data/pigeon.fasta"),
    "rabbit": Sequence("data/rabbit.fasta")
}

print("\n\n")
print("Creating Phylogenic Tree:")
print("===============================================================")
print("Input Sizes:")
for s in species:
    print("Length of sequence for species", s, "is:", species[s].getSize())

# Printing Results
tree = create_phylogenic_tree(species)
print("\nTree Representation ([A, B] means A and B are linked genetically)")
print("\nOutput:")
print(tree)
print("===============================================================\n")

"""
    [[elephant, [frog, dog]], [rabbit, [pigeon, [human, whale]]]]
    means 
           __human
        __| 
     __|  |__whale
  __|  |__pigeon
 |  |__rabbit
 |
 |      __frog
 |   __|
 |__|  |__dog
    |__elephant

"""