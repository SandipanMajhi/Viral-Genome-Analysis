import argparse
from random import shuffle
from tqdm import tqdm 
import json

class Similarity:
    def __init__(self, file1, file2, k):
        '''
        file_list : list of the fasta files
        '''
        self.file1 = file1
        self.file2 = file2
        self.kmers = {}
        self.union_set = set()
        self.k = k
        self.read1 = ""
        self.read2 = ""
        self.read_file()
        self.kmer1 = self.build_kmers(self.read1)
        self.kmer2 = self.build_kmers(self.read2)
        self.union_set = self.build_union()
        self.inter_set = self.build_intersection()
        self.jacard = self.jacard_index()


    def read_file(self):
        with open(self.file1, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                if line[0] != ">":
                    self.read1 += line
        
        with open(self.file2, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                if line[0] != ">":
                    self.read2 += line
        

    def build_kmers(self, read):
        # print("Building Table")
        # print("Read length ------ {} and {} ------".format(len(self.read1), len(self.read2)))
        i = 0
        j = i + self.k
        rsize = len(read)
        kmer_set = set()
        while j <= rsize:
            kmer = read[i:j]
            kmer_set.add(kmer)
            i += 1
            j = i + self.k
        return kmer_set

    def build_union(self):
        return self.kmer1.union(self.kmer2)
        
    def build_intersection(self):
        return self.kmer1.intersection(self.kmer2)

    def jacard_index(self):
        return len(self.inter_set)/len(self.union_set)

    def overlap_coefficient(self):
        return len(self.inter_set)/min(len(self.kmer1), len(self.kmer2))

def main():
    ### First get the kmers
    ap = argparse.ArgumentParser()
    ap.add_argument("-file1", default="data/inA_Sh_H7N9.fna", type=str, help="input first viral genome")
    ap.add_argument("-file2", default="data/inA_HK_H9N2.fna", type=str, help="input second viral genome")
    ap.add_argument("--k", default=32, type = int, help = "input k-mer size")
    ap.add_argument("-o", default="output", type = str, help= "output folder location")
    ap.add_argument("--n", default=150, type=int, help= "k-mer limit")
    args = ap.parse_args()

    for ks in range(args.k, args.n+1 ):
        J = Similarity(args.file1, args.file2, ks)
        print(J.jacard)


if __name__ == "__main__":
    main()