# Viral-Genome-Analysis

### Please donot run any other py files.

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
