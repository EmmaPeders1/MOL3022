﻿# MOL3022 - SVE Protein Alignment Tool

The SVE (Solveig, Vebjørn and Emma) Protein Alignment Tool is a tool that finds the best protein alignments for a user-given query sequence.
The project is made for the course [MOL3022 - Bioinformatics - Method Oriented Project](https://www.ntnu.edu/studies/courses/MOL3022#tab=omEmnet) and consists of a Python script with a Streamlit-based front-end.

## How to run the project

1. Download protein sequences in FASTA (canonical) format. We recommend downloading protein sequences from [UniProtKnowledgeBase](https://www.uniprot.org/uniprotkb?query=reviewed:true) by clicking the "Download" button. You may either download all files or choose some yourself. If you download the files as compressed make sure you unzip the files after download.

<img src="docs/Docs_UniProtKBDownloadHighlighted.png" width="800">

2. Create a folder named `uniprot_fasta` in this repository and put your fasta file in it. Make sure the fasta file is named `uniprot_sprot.fasta` which is the default naming when downloading from UniProtKB.

3. Open a terminal and use the command `pip install -r requirements.txt`. Make sure your terminal is in the `MOL3022` folder.

4. Run the project using `streamlit run main.py` in your terminal.

## Expected results

When running the project using Streamlit your browser should open and connect to localhost:8051. (or a different port if 8501 is in use)
After loading you will be met by the following page:

<img src="docs/SVE_Main_Page.png" width="800">

Here you can choose how many sequences you want to compare alignment to. This is useful when you have a lot of sequences (such as if you download all sequences from UniProtKB) as comparing alignment to all of them can take a lot of time and may not be necessary when just testing the software. The software also limits you to one best alignment from each sequence. This is because one sequence can give huge amounts of possible best alignments and this scales exponentially with the amount of alignments you compare to. By giving one alignment for each sequence, a user can instead look from there at which sequence they would like to study further.

After selecting the amount of sequences you want to compare alignments to you can choose to input a query sequence or use the default that is already input and then click the `Find closest alignment` button.

The software will then begin analyzing your query sequence to each other sequence and the progress bar indicates how many sequences it has completed alignment to compared to how many are left.
<img src="docs/Docs_ProgressBar.png">

After the best alignments (from unique sequences) have been found they are presented as a list.
<img src="docs/Alignments_Found.png" width="800">

## Scoring

Our tool uses BLOSUM62 as our scoring matrix. BLOSUM is widely regarded as a good scoring matrix for detecting most weak protein similarties.
