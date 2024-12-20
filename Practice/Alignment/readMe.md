# Seqeunce Alignment and Phylogeny

Sequence Alignment is crucial bioinformatics task and one of the tools for this task is [Clustal Omega](http://www.clustal.org/omega/)

### Terminal


If you are new to command line and terminal please watch [this video](https://www.youtube.com/watch?v=jgcXclSXnVo&t=44s)

### Installing CLustalO


Now open the terminal and let's start:

To install in your machine:
```
sudo apt install clustalo
```

To test the installation:
```
clustalo -h
```


Create a directory in home named as 'Projects'
```
cd ~
```
Now you are in home
```
mkdir Projects
```

```
cd Projects
```
Now you are in Projects

### Fasta File Format

There are very well defined file formats in the bioinformatics to represent the data, and knowing them is crucial. One of them is .fasta format. Please visit [this site](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) 

Create and empty file:
```
touch Pseudomonas_VD2_16S_rRNA.fasta
```

And open it, add this 16S sequence:

```
>Pseudomonas_VD2_16S_rRNA_gene
CTAGACAGTGGCTATTTCATCAGCCGTCAGTTCACGAAACTCCCCGGGCGCCAAACCAGGGTCCAGACAGATTGCGCCCATGCTTTCCCGGTGCAACCCCACCACCTTGTTGTTGAAGTGCCCGAACATGCGCTTGACCTGGTGATAACGCCCTTCGACGATCGCCAGCCGGGCCTGGCGCGGGCCGAGGATGTCCAGCAGGGCAGGTTGGGTGGTGAGGTTTTCGAAGGCGAAATAGAACCCCTCGCGAAACCTGGCCACATAGTGCTCGCCAATTTCGTCCTCGGTGTCCACCAGGTAATGCTTGGGCAACTTGGTGGCAGGCTGGGTCAGGCGCCTTGACCACTGGCCATCGTTGGTCAGGATCATCAGGCCGGTGGTGTTGAAGTCCAGGCGCCCGGCTATGTGCAGGTCATCCCGCAACGCTGCGGGCAGCAGGTCGAGCACGGTCGGGTGCTGTGGGTCTTGGGTGGCGCTGACGCAGCCAGTGGGCTTGTGCAGCATCAGGTAGCGGGCCGGGCGGCCGGCCTGCAGCAGTTGCTCGTCCAGCTCCACGCGGCTGAATTCACGCACTTCGGCCAGTGGGTCGCTGACAACCTGGCCGTCCACCCGCACGCGGCGTTGCGCCAGCATCAGTCGGGCTTGCTGGCGGTTGTGGCTGGGCAGGTTGGCGAGAAAACGGTCAAGGCGCAT```
```

### Blast
Now lets find other closest species;

Blast is very powerful search engine (some says it is Google of biology), 

[For info:](https://youtu.be/WRKQGwh_Mw0?si=wVb43f8I8gB8FWQr)

Open the link:
https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastn&PAGE_TYPE=BlastSearch&LINK_LOC=blasthome

And search your 16S_rRNA sequence with blastn.

in the results go to **Alignments** tab, download first 10 results with the option **Download**, **FASTA, aligned sequence**

Now, rename the downloaded files, and combine them into one fasta file named as **Combined_Ecoli.fasta**:

### Alignment
Before going further, What is [alignment](https://www.youtube.com/watch?v=LhpGz5--isw)? 


And run the clustalo alignment:
```
clustalo -i Combined_Ecoli.fasta -o Pseudomonas_VD2_16S_rRNA.aln --outfmt=clu --guidetree-out=tree_Ecoli_16S_rRNA.dnd
```

Now, you can see the results are in different formats, search this format types in the google and find some way to visualize .dnd format. What is .aln format?


### Hints:
  - Carl Richard Woese
  - https://itol.embl.de/upload.cgi
  - https://www.youtube.com/watch?v=KTZIYt8z9zQ


### Question:
  - What is the **.sp** meaning in the naming of species: Pseudomonas sp. VD2 
  - What is the species level classification of Pseudomonas sp. VD2
  - Why we use 16S rRNA for species determination?
  - Why the sequence alignment is important?
  - What Blast do?





