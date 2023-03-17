import argparse
from Kmer import Kmer
from sim import Sim
from tqdm import tqdm
import json

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-file1", default="data/inA_Sh_H7N9.fna", type=str, help="input first viral genome")
    ap.add_argument("-file2", default="data/inA_HK_H9N2.fna", type=str, help="input second viral genome")
    ap.add_argument("--k", default=32, type = int, help = "input k-mer size")
    ap.add_argument("--m", default=4, type = int, help = "Allowed Mismatch")
    ap.add_argument("-o", default="output", type = str, help= "output folder location")
    ap.add_argument("--n", default=50, type=int, help= "k-mer limit")
    args = ap.parse_args()
    outfile1 = args.file1.split("/")[1].split(".")[0] + ".json"
    outfile2 = args.file2.split("/")[1].split(".")[0] + ".json"
    match_file = "kmers_HK_Sh.json"
    
    kmer_ranges = {}

    for ks in tqdm(range(args.k, args.n+1)):
        kmers = Kmer(args, ks)
        sim = Sim(outfile1, outfile2)
        if sim.num_matches != 0:
            kmer_ranges[ks] = sim.num_matches
        else:
            break
        # merge = Merge_Kmer(match_file, args.k, kmers.read1, kmers.read2)
    
    with open(match_file, "w") as fp:
        json.dump(kmer_ranges, fp)


if __name__ == "__main__":
    main()

