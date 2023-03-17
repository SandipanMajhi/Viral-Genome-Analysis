import json
import itertools
import shutil

class Sim:
    def __init__(self, file1, file2):
        self.folder = "./output/"
        self.f1 = file1
        self.f2 = file2
        self.d1 = None
        self.d2 = None 
        self.num_matches = 0
        self.read_files()
        self.find_matches()
        
    def read_files(self):
        with open(self.folder+self.f1, "r") as fp:
            self.d1 = json.load(fp)
        with open(self.folder+self.f2, "r") as fp:
            self.d2 = json.load(fp)

    def find_matches(self):
        self.matches = {}
        for k,v in self.d1.items():
            if k in self.d2:
                list_product = itertools.product(self.d1[k], self.d2[k])
                self.matches[k] = []
                for elem in list_product:
                    self.matches[k].append(elem)
        self.num_matches = len(self.matches)
        file_name = "match.json"
        with open(file_name, "w") as fp:
            json.dump(self.matches, fp, indent=1)
        
        src = file_name
        dest = "./output/" + src
        shutil.move(src, dest)