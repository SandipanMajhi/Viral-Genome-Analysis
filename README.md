# Viral-Genome-Analysis

### Main Aim ###

In this project we take multiple strains of Influenza A virus and trace how they spread over-time. Usually viral strains travel through people and we trace this phenomena using vector based techniques and show how this can be very useful.

### Total Explanation with PDF ### 

Go to this link to understand the implementation of the algorithm https://drive.google.com/drive/folders/1FRIowBfnnT5RB79RBAZ8WdPv_29J6VY7?usp=sharing 

### Results of the Phylogenetic Experiment ###
![s1](https://github.com/SandipanMajhi/Viral-Genome-Analysis/assets/52757231/8e4f3756-b354-4056-9857-4d7d03963a5f)

![s2](https://github.com/SandipanMajhi/Viral-Genome-Analysis/assets/52757231/1ad4fe7a-a2ea-4106-8da0-436bea3d413e)




### Please donot run any other py files except mentioned in the following.


main.py  - generates all the kmer frequencies with two genomes

    -file1 : input first viral genome
    -file2 : input second viral genome
    --k : input k-mer size
    --m : allowed mismatch
    -o : output file location
    --n : kmer length limit
    

jacard_full.py : calculates the jacard coefficients over all the sequences in a directory

    -path : path to the directory of the genomes
    --k : kmer search starts
    --max_k : max kmer size to be searched
    -mm_prob : allowed mismatch prob
    --n : number of kmers to be sampled for vector representation
    

cosine_vec.py : calculates the cosine similairty over all the sequences in a directory

    -path : path to the directory of the genomes
    --k : kmer search starts
    --max_k : max kmer size to be searched
    -mm_prob : allowed mismatch prob
    --n : number of kmers to be sampled for vector representation
