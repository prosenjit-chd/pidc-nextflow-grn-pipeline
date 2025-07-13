import matplotlib.pyplot as plt
import scanpy as sc

# Load input genes from h5ad
adata = sc.read_h5ad("filtered_placeholder.h5ad")
if adata is None or adata.var_names is None:
    raise ValueError("❌ 'adata.var_names' is missing or empty.")

input_genes = set(adata.var_names)

# Load output genes from the CSV
output_genes = set()
with open("results/pidc_output_final_rankedEdges.csv", "r") as file:
    for line in file:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            output_genes.add(parts[0].strip())
            output_genes.add(parts[1].strip())

# Compare
only_input = input_genes - output_genes
only_output = output_genes - input_genes
common_genes = input_genes & output_genes

# Plot bar chart
plt.figure(figsize=(8, 5))
plt.bar(
    ["Only in Input", "Only in Output", "Common Genes"],
    [len(only_input), len(only_output), len(common_genes)],
    color=["skyblue", "salmon", "lightgreen"],
)
plt.ylabel("Gene Count")
plt.title("Gene Comparison: Input (.h5ad) vs Output (.csv)")
plt.grid(axis="y")
plt.tight_layout()
plt.savefig("gene_comparison.png")
plt.show()
