#### Similarity measure over a range of ks

import argparse
import json
import random
from random import shuffle
from tqdm import tqdm
import os
import numpy as np
from numpy.linalg import norm
import itertools


class Similarity:
    def __init__(self, file_list, k, max_k, mm_prob, num_samples):
        self.file_list = file_list
        self.num_samples = num_samples
        self.mm_prob = mm_prob
        self.sequence_vectors = {}
        self.k = k
        self.max_k = max_k
        self.similarity_scores = {}
        self.kmers = {}
        self.reads = {}
        self.get_reads()
        self.get_kmers()
        self.union_set = set()
        self.build_union_set()
        self.build_vectors(self.mm_prob)
        self.get_jacard_similarity()
        self.get_relation_matrix()

    def get_reads(self):
        print("Getting the reads.....")
        for file in self.file_list:
            self.kmers[file] = set()
            read = ""
            with open(file, "r") as fp:
                for line in fp.readlines():
                    line = line.strip()
                    if line[0] != ">":
                        read += line
            self.reads[file] = read

    def get_kmers(self):
        print("Forming kmers.....")
        for file in self.file_list:
            if file not in self.kmers:
                self.kmers[file] = set()
        
        for file in self.file_list:
            for ks in tqdm(range(self.k, self.max_k)):
                i = 0
                j = i + ks
                read_len = len(self.reads[file])
                while j <= read_len:
                    self.kmers[file].add(self.reads[file][i:j])
                    i += 1
                    j = i + ks
            
    def build_union_set(self):
        print("Forming union set.....")
        for file in self.file_list:
            self.union_set = self.union_set.union(self.kmers[file])

    def build_bi_union_set(self, pairs):
        # print("Forming the Bi Union set.....")
        bi_union_set = self.kmers[pairs[0]].union(self.kmers[pairs[1]])
        return bi_union_set


    def get_hamming_distance(self, s1, s2):
        dist = 0
        for i in range(len(s1)):
            if s1[i] != s2[i]:
                dist += 1
        return dist

    def build_intersections(self):
        all_pairs = list(itertools.combinations(self.file_list, 2))
        self.intersects = {}
        for i in tqdm(range(len(all_pairs))):
            self.intersects[all_pairs[i]] = self.kmers[all_pairs[i][0]].intersection(self.kmers[all_pairs[i][1]])
            lf1 = list(self.kmers[all_pairs[i][0]])
            lf2 = list(self.kmers[all_pairs[i][1]])
            for km1 in lf1:
                for km2 in lf2:
                    if len(km1) == len(km2):
                        ham_dist = self.get_hamming_distance(km1, km2)
                        if ham_dist/len(km1) <= self.mm_prob:
                            self.intersects[all_pairs[i]].add(km1)
                            self.intersects[all_pairs[i]].add(km2)


    def build_vectors(self, mm_max = 0.2):
        print("Building vectors.....")
        union_list = list(self.union_set)
        shuffle(union_list)
        sample_kmers = random.sample(union_list, self.num_samples)
        self.union_list = sample_kmers
        vectors = {}

        for kmer in union_list:
            vectors[kmer] = {}

        for file in self.file_list:
            for kmer in tqdm(self.union_list):
                vectors[kmer][file] = 0
                if kmer in self.kmers[file]:
                    vectors[kmer][file] = 1
                else:
                    for km in list(self.kmers[file]):
                        if len(km) == len(kmer):
                            ham_dist = self.get_hamming_distance(km, kmer)
                            if ham_dist/len(kmer) <= mm_max:
                                vectors[kmer][file] = 1 

        for file in self.file_list:
            self.sequence_vectors[file] = np.zeros((len(self.union_list), 1))
            for i in range(len(self.union_list)):
                self.sequence_vectors[file][i,0] = vectors[self.union_list[i]][file]

    def get_jacard_similarity(self):
        all_pairs = itertools.combinations(self.file_list, 2)
        for pairs in all_pairs:
            self.similarity_scores[pairs] = np.sum(np.multiply(self.sequence_vectors[pairs[0]], self.sequence_vectors[pairs[1]])) / len(self.union_list)
    
    def get_relation_matrix(self):
        file_dict = {}
        i = 0
        for file in self.file_list:
            file_dict[file] = i
            i += 1

        all_pairs = itertools.combinations(self.file_list, 2)
        self.sim_mat = np.zeros((len(self.file_list), len(self.file_list)))
        for pairs in all_pairs:
            self.sim_mat[file_dict[pairs[0]] , file_dict[pairs[1]]] = self.similarity_scores[pairs]




def get_all_files(path):
    file_list = []
    for (dirpath, dirname, fnames) in os.walk(path):
        for file in fnames:
            if file.startswith("in"):
                file_list.append(dirpath+os.sep+file)
    return file_list


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-path", default="influenza", type=str, help="data path")
    ap.add_argument("--k", default=75, type = int, help = "start k mer size search")
    ap.add_argument("--max_k", default=150, type = int, help = "max kmer size")
    ap.add_argument("-mm_prob", default=0.1, type = float, help= "mismatch rate")
    ap.add_argument("--n", default=100, type=int, help="num samples of kmers")

    args = ap.parse_args()
    genome_files = get_all_files(args.path)
    Sim = Similarity(genome_files, k=args.k, max_k=args.max_k, mm_prob = args.mm_prob, num_samples=args.n)
    print(Sim.similarity_scores)

if __name__ == "__main__":
    main()
