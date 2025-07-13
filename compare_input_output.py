import scanpy as sc
import pandas as pd

# Load input genes from h5ad
adata = sc.read_h5ad("filtered_placeholder.h5ad")
input_genes = set(adata.var_names)

# Load output gene interactions from PIDC
edges_df = pd.read_csv(
    "results/pidc_output_final_rankedEdges.csv", sep="\t", header=None
)
output_genes = set(edges_df[0].unique()).union(set(edges_df[1].unique()))

# Compare
common_genes = input_genes.intersection(output_genes)
only_in_input = input_genes - output_genes
only_in_output = output_genes - input_genes

print(f"✅ Input gene count: {len(input_genes)}")
print(f"✅ Output gene count: {len(output_genes)}")
print(f"✅ Common genes: {len(common_genes)}")
print(f"❌ Genes only in input: {len(only_in_input)}")
print(f"❌ Genes only in output: {len(only_in_output)}")

# Optionally save result
with open("comparison_summary.txt", "w") as f:
    f.write(f"Input genes: {len(input_genes)}\n")
    f.write(f"Output genes: {len(output_genes)}\n")
    f.write(f"Common genes: {len(common_genes)}\n")
    f.write(f"Only in input: {only_in_input}\n")
    f.write(f"Only in output: {only_in_output}\n")
