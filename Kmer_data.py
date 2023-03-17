import json
import shutil

class Kmer:
    def __init__(self, args):
        self.f1 = args.file1
        self.f2 = args.file2
        self.k = args.k
        self.output_folder = args.o
        self.read1 = ""
        self.read2 = ""
        self.read_file()
        self.store_dict(self.build_table(self.read1), self.f1.split("/")[-1].split(".")[0])
        self.store_dict(self.build_table(self.read2), self.f2.split("/")[-1].split(".")[0])

    def read_file(self):
        print("Reading files --------------{} and {}".format(self.f1, self.f2))
        with open(self.f1, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                if line[0] != ">":
                    self.read1 += line
        
        with open(self.f2, "r") as fp:
            for line in fp.readlines():
                line = line.strip()
                if line[0] != ">":
                    self.read2 += line

    def build_table(self, read):
        print("Building Table")
        print("Read length ------ {} and {} ------".format(len(self.read1), len(self.read2)))
        i = 0
        j = i + self.k
        rsize = len(read)
        kmer_dict = {}

        while j <= rsize:
            kmer = read[i:j]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = [i]
            else:
                kmer_dict[kmer].append(i)
            i += 1
            j = i + self.k
        return kmer_dict

    def get_closest_unique_table(self, read):
        ks = self.k
        pass


    def store_dict(self, kmer_dict, file_name):
        print("Storing Dicts -----------------")
        with open(file_name+".json", "w") as fp:
            json.dump(kmer_dict, fp)

        src = file_name + ".json"
        dest = "./output/" + src
        shutil.move(src, dest)
