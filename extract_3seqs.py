seqs_to_keep = [
    "gi|568815592:31353871-31357211",  # Primary 
    "gi|299782605:5000-8340",           # 100% identity (RC)
    "gi|568815529:2834231-2837570"      # 98.7% identity
]

keep = False
for line in open("HLA-zoo/seqs/B-3106.fa"):
    if line.startswith(">"):
        seq_name = line.split()[0][1:]
        keep = seq_name in seqs_to_keep
    if keep:
        print(line, end='')
