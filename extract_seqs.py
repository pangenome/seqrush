import sys

seqs_to_keep = [
    "gi|568815592:31353871-31357211",
    "gi|299782605:5000-8340"
]

keep = False
for line in open("HLA-zoo/seqs/B-3106.fa"):
    if line.startswith(">"):
        seq_name = line.split()[0][1:]
        keep = seq_name in seqs_to_keep
    if keep:
        print(line, end='')
