```python
!pip install biopython
from Bio import SeqIO

#read in fasta file and convert it to string
rna = ""
from Bio import SeqIO
for seq_rec in SeqIO.parse("sequence.fa", "fasta"):
    seq_str = str(seq_rec.seq)
for i in seq_str:
    if i == "A":
        rna += "A"
    if i == "C":
        rna += "C"
    if i == "G":
        rna += "G"
    if i == "T":
        rna += "U"
print(rna[:100])

#compute reverse complement
from Bio import SeqIO
for seq_rec in SeqIO.parse("sequence.fa", "fasta"):
  rc = seq_rec.seq.reverse_complement()
  print(str(rc)[:100])

#build dictionary of k-mer frequencies for k=3 and k=4
from Bio import SeqIO
for seq_rec in SeqIO.parse("sequence.fa", "fasta"):
    seq_str = str(seq_rec.seq)
threemers = {}
for i in range(len(seq_str)-2):
    threemer = seq_str[i:i+3]
    if threemer in threemers:
        threemers[threemer] += 1
    else:
        threemers[threemer] = 1
fourmers = {}
for i in range(len(seq_str)-3):
    fourmer = seq_str[i:i+4]
    if fourmer in fourmers:
        fourmers[fourmer] += 1
    else:
        fourmers[fourmer] = 1
  
#plot top 50 most frequent 4-mers
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
df = pd.DataFrame(list(fourmers.items()), columns=["4-mer", "Frequency"])
df = df.sort_values(by="Frequency", ascending=False)
top_50_df = df.head(50)
top_50_df["4-mer"] = pd.Categorical(top_50_df["4-mer"], categories=top_50_df["4-mer"], ordered=True)
plt.figure(figsize=(12,8))
sns.barplot(data=top_50_df, x="4-mer", y="Frequency")
plt.title("Top 50 Most Frequent 4-mers", fontsize=16)
plt.xlabel("4-mers", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.xticks(rotation=90, fontsize=10)
plt.tight_layout()
plt.show()

#plot all 3-mers
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
df = pd.DataFrame(list(threemers.items()), columns=["3-mer", "Frequency"])
df = df.sort_values(by="Frequency", ascending=False)
threemers = df.head(64)
threemers_top64["3-mer"] = pd.Categorical(threemers_top64["3-mer"], categories=threemers_top64["3-mer"], ordered=True)
plt.figure(figsize=(12,8))
sns.barplot(data=threemers, x="3-mer", y="Frequency")
plt.title("Frequency of 3-mers", fontsize=16)
plt.xlabel("3-mers", fontsize=14)
plt.ylabel("Frequency", fontsize=14)
plt.xticks(rotation=90, fontsize=10)
plt.tight_layout()
plt.show()

#convert mRNA sequences from FASTA file into proteins
from Bio import SeqIO
for line in SeqIO.parse("transcriptome.fa", "fasta"):
    if line.id == "lcl|NC_045512.2_gene_2":
        mrna = line.seq
        s_protein = mrna.translate(to_stop=False)
print(s_protein)

#Create normalized amino acid composition bar chart per transcript
from Bio import SeqIO
import pandas as pd
import matplotlib.pyplot as plt

data = []
for seq in SeqIO.parse("transcriptome.fa", "fasta"):
    gene = seq.id
    trim = seq.seq[:len(seq.seq) - len(seq.seq) % 3]
    aa_seq = trim.translate(to_stop=False)
    for position, aa in enumerate(aa_seq, start=1):
        data.append({"Gene": gene, "Position": position, "AA": aa})
df = pd.DataFrame(data)
pivot_df = df.pivot_table(index="Gene", columns="AA", aggfunc="size", fill_value=0)
normalized_df = pivot_df.div(pivot_df.sum(axis=1), axis=0)
ax = normalized_df.plot(kind="bar", stacked=True, figsize=(12, 8), colormap="tab20")
plt.title("Normalized Amino Acid Composition Per Transcript", fontsize=16)
plt.xlabel("Transcript", fontsize=14)
plt.ylabel("Proportion", fontsize=14)
plt.xticks(rotation=90, fontsize=10)
plt.legend(title="Amino Acids", bbox_to_anchor=(1.05, 1), loc="upper left")
plt.tight_layout()
plt.show()
```
