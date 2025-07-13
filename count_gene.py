import pandas as pd

# Use correct path to the file
df = pd.read_csv("results/pidc_output_final_rankedEdges.csv", sep="\t", header=None)
genes = set(df[0]).union(set(df[1]))
print(f"Total unique genes: {len(genes)}")
