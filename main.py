from Bio import SeqIO
from Bio import Align
import streamlit as st
import re
import numpy as np

# SVE (Solveig, Vebjørn and Emma) Alignment Tool
# This tool is designed to find the closest alignment to a query sequence in a fasta file
# Use pip install -r requirements.txt to install required packages

# Maximum number of sequences to align
sequence_cap = 10000

# Maximum number of alignments for each sequence
alignment_cap = 50000

# Regex pattern to find sequence name in description
pattern = r'OS=(.*?)\sOX='

# Set page config
st.set_page_config(layout="wide")

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
    # Find the closest alignment to the query sequence in the sequences

    # Stores aligment parameters
    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    best_alignments = []
    best_score = float("-inf")
    print("Finding best alignment, please wait...")

    # Find alignment with highest score
    counter = 0

    # Progress bar
    progress_text = "Finding best alignment. Please wait..."
    progress_bar = st.progress(0, text=progress_text)

    for identifier, sequence in sequences.items():
       
        try:
            # Use sequence[0] as sequence[1] is the name
            alignments = aligner.align(sequence[0], query_sequence) #Alignment objects represtenting alignments between sequence and query.
        except: # If the alignment score is too high, it will cause an overflow error, so we throw it away and continue
            continue
            
        # Create for-loop counter and show progress
        counter += 1
        progress_bar.progress(counter / len(sequences), text=progress_text)
        
        try:
            
            if 0 < len(alignments) < alignment_cap:
                best_alignment = None                
                for alignment in alignments:
                    if alignment.score > best_score:
                        best_score = alignment.score
                        best_alignment = alignment
                        
                if best_alignment:
                    best_alignments.append(best_alignment)
        except: # If the length of the alignments is too high, it will cause an overflow error, but that also indicates that the alignment is not useful, so we throw it away and continue
            continue
    
    # Remove progress bar and show completion
    progress_bar.empty()
    st.balloons()

    st.markdown(":blue[Alignment complete!]")
    result = [best_alignments, best_score]
    return result

if __name__ == "__main__":
    st.header("SVE Protein Alignment Tool :rocket:")
    st.markdown("*This is an alignment tool that finds the closest protein alignments in a fasta file to a given query sequence*")
    fasta_file = "uniprot_fasta/uniprot_sprot.fasta"
    st.markdown(":blue[Reading fasta file: " + fasta_file + "]")
    sequences = read_fasta(fasta_file)
    seq_amount = len(sequences)
    st.markdown("**Number of Sequences:** " + str(seq_amount))

    # Spacing 
    st.text("")

    # Give user option to use a subset of sequences if there are more than 50k sequences
    if seq_amount > sequence_cap:
        st.markdown(f":red[WARNING: More than {sequence_cap} sequences! Recommended to use a subset.]") # With all 500k sequences, the program can a long time to run
        # Give user option to use their chosen amount of sequences
        seq_slider = st.slider("Select the wanted amount of sequences: ", 0, seq_amount, sequence_cap, 50)
        sequences = {k: sequences[k] for k in list(sequences)[:seq_slider]}
    st.markdown("**" + str(seq_slider) + " sequences will be used for alignment.**")    

    # Spacing 
    st.text("")

    # Give user option to enter query sequence
    user_sequence = st.text_input("Enter query sequence (or use the default sequence in place): ", "MAFSAEDVLKEYDRRRRMEALLLSLYYPNDRKLLDYKEWSPPRVQVECPKAPVEWNNPPEKGLIVGHFSGIKYKGEKAQASEVDVNKMCCWVSKFKDAMRRYQGIQTCKIPGKVLSDLDAKIKAYNLTVEGVEGFVRYSRVTKQHVAAFLKELRHSKQYENVNLIHYILTDKRVDIQHLEKDLVKDFKALVESAHRMRQGHMINVKYILYQLLKKHGHGPDGPDILTVKTGSKGVLYDDSFRKIYTDLGWKFTP")
    
    # Start alignment
    if st.button("Find Closest Alignment", type="primary"):
        resultList = find_closest_alignment(sequences, user_sequence)
        len_result = len(resultList[0])
        st.markdown(f"**{len_result} alignment(s) found with score: {str(resultList[1])}**")
        
        # Show results
        st.subheader("Best Alignment(s):")
        for i in range(len(resultList[0])):
            st.markdown(f"{sequences[list(sequences.keys())[i]][1]} - {list(sequences.keys())[i]}")
            st.code(str(resultList[0][i]))