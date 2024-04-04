from Bio import SeqIO
from Bio import pairwise2
from Bio import Align

# SVE (Solveig, Vebjørn and Emma) Alignment Tool
# This tool is designed to find the closest alignment to a query sequence in a fasta file
# Use pip install -r requirements.txt to install required packages

def read_fasta(filename):
    # Read fasta file and return dict of sequences
    sequences = {}
    for record in SeqIO.parse(filename, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def find_closest_alignment(sequences, query_sequence):
    # Find closest alingment to query sequence using globalxx alignment
    # TODO: find alignment witg highest score takes ages, so needs to be fixed (broken)
    
    #stores aligment parameters
    aligner = Align.PairwiseAligner()

    best_alignment = []
    best_score = float("-inf")
    print("Finding best alignment, please wait...")

    #find highest score
    for identifier, sequence in sequences.items():

        #calculte the alignment score
        score = aligner.score(sequence, query_sequence)


        if score > best_score:
            best_score = score

    #find alignment with highest score
    for identifier, sequence in sequences.items():    
        alignments = aligner.align(sequence, query_sequence) #Alignment objects represtenting alignments between sequence and query

        for x in range(len(alignments)):
            if alignments[x].score == best_score:
                best_alignment.append(alignments[x].score)
    

       
    result = [len(best_alignment), best_score]
    return result
# Main
print("\n SVE Alignment Tool\n---------------------\n")
fasta_file = "uniprot_fasta/uniprot_sprot.fasta"
print("Reading fasta file:", fasta_file, "please wait...")
sequences = read_fasta(fasta_file)
print("Number of Sequences:", len(sequences), "\n")
if len(sequences) > 50000:
    print('\x1b[1;31;40m' + "WARNING: More than 50k sequences! Recommended to use a subset for testing." + '\x1b[0m') # With all 500k sequences, the program can take several minutes to run
    if input("Use first 50 sequences? (y/n): ").lower() == "y":
        sequences = {k: sequences[k] for k in list(sequences)[:50]}
print(len(sequences), "sequences will be used for alignment.")
query_sequence = "REHSYWDSWSHKSMWYDDGCACPFGNNLHFHHPWANNYSCLTRIKFVIFM" # Example query sequence, can be replaced with input from user
resultList = find_closest_alignment(sequences, query_sequence)


print("Best Alignment:")
#print("Identifier:", best_alignment[0])
print("Sequence:")
print(resultList[0])

print("Alignment Score:", resultList[1])
#print("Alignment:")
#print(pairwise2.format_alignment(*best_alignment[2]))
