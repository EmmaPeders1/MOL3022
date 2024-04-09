from Bio import SeqIO
from Bio import pairwise2
from Bio import Align
import re

# SVE (Solveig, Vebjørn and Emma) Alignment Tool
# This tool is designed to find the closest alignment to a query sequence in a fasta file
# Use pip install -r requirements.txt to install required packages

sequence_cap = 10000
pattern = r'OS=(.*?)\sOX='

def read_fasta(filename):
    # Read fasta file and return dict of sequences
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        match = re.search(pattern, record.description)
        if match:
            sequences[record.id] = [str(record.seq), match.group(1)]
        else:
            sequences[record.id] = [str(record.seq), record.description]
        
    return sequences

def find_closest_alignment(sequences, query_sequence):
    # Find closest alingment to query sequence using globalxx alignment
    # TODO: find alignment witg highest score takes ages, so needs to be fixed (broken)
    
    #stores aligment parameters
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 2
    aligner.mismatch_score = -2
    aligner.open_gap_score = -0.5
    aligner.extend_gap_score = -0.1
    

    best_alignment = []
    best_score = float("-inf")
    print("Finding best alignment, please wait...")

    

    #find alignment with highest score
    counter = 0
    for identifier, sequence in sequences.items():
        try:
            alignments = aligner.align(sequence[0], query_sequence) #Alignment objects represtenting alignments between sequence and query
        except OverflowError: # Skip if alignment is too long
            continue
        # Create for loop counter and print
        counter += 1
        print("Aligning sequence", counter, "of", len(sequences), end="\r")

        if 0 < len(alignments) < 50000: 
            for x in range(len(alignments)):
                if alignments[x].score >= best_score:
                    # print(alignments[x])
                    if alignments[x].score > best_score:
                        best_score = alignments[x].score
                        best_alignment = []
                    best_alignment.append(alignments[x])
                    print(alignments[x])
                    print(sequence[1])
                    # Only keep one best alignment from each sequence
                    break
        
    print("\nAlignment complete.")
       
    result = [best_alignment, best_score]
    return result
# Main
print("\n SVE Alignment Tool\n---------------------\n")
fasta_file = "uniprot_fasta/uniprot_sprot.fasta"
print("Reading fasta file:", fasta_file, "please wait...")
sequences = read_fasta(fasta_file)
print("Number of Sequences:", len(sequences), "\n")
if len(sequences) > sequence_cap:
    print('\x1b[1;31;40m' + "WARNING: More than " + str(sequence_cap) + " sequences! Recommended to use a subset for testing." + '\x1b[0m') # With all 500k sequences, the program can take several minutes to run
    if input("Use first " + str(sequence_cap) + " sequences? (y/n): ").lower() == "y":
        sequences = {k: sequences[k] for k in list(sequences)[:sequence_cap]}
print(len(sequences), "sequences will be used for alignment.")
query_sequence = "REHSYWDSWSHKSMWYDDGCACPFGNNLHFHHPWANNYSCLTRIKFVIFM" # Example query sequence, can be replaced with input from user
resultList = find_closest_alignment(sequences, query_sequence)


# Print results
print("Best alignment: ")
print(resultList[0][0])
print("______________________________________")
print("Number of alignments with best score: ")
print(len(resultList[0]))
print("______________________________________")
print("Best alignment sequence ")
print(resultList[0][0][0])

print("Alignment Score:", resultList[1])


