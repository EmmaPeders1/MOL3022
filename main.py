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
    # TODO: Choose correct alignment in list of alignment objects
    
    #stores aligment parameters
    aligner = Align.PairwiseAligner()

    best_alignment = None
    best_score = float("-inf")
    print("Finding best alignment, please wait...")

    for identifier, sequence in sequences.items():
       #alignments = pairwise2.align.globalxx(query_sequence, sequence, one_alignment_only=True)
        #alignment = alignments[0]
        #score = alignment[2]

        #calculte the alignment score
        score = aligner.score(sequence, query_sequence)


        if score > best_score:
            best_score = score
            best_alignment = aligner.align(sequence, query_sequence)[0] #Alignment objects represtenting alignments between sequens and query

       # if score > best_score:
           # best_score = score
           # best_alignment = (identifier, sequence, alignment)
    result = [best_alignment, best_score]

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
for resultAlign in resultList[0]:
    print("Sequence:")
    print(resultAlign)

print("Alignment Score:", resultList[1])
#print("Alignment:")
#print(pairwise2.format_alignment(*best_alignment[2]))
