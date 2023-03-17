import json
import matplotlib.pyplot as plt

save_dir = "plot/"

file_path = "kmers_UK_Gu.json"

with open(file_path, "r") as fp:
    data = json.load(fp)

x = []
y = []

for k,v in data.items():
    x.append(k)
    y.append(v)

plt.plot(x,y)
plt.title("Matching Kmers between UK-H1N1(2010) and Gu-H5N1(1996)")
plt.xticks(x[::10])
plt.xlabel("K-MER SIZE")
plt.ylabel("MATCH FREQUENCY")
plt.savefig("UK_Gu.png")
plt.show()