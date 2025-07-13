import matplotlib.pyplot as plt
import networkx as nx

input_file = "results/pidc_output_final_rankedEdges.csv"
edges = []

# Manually read and parse the file line-by-line
with open(input_file, "r") as file:
    for line in file:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            try:
                gene1 = parts[0].strip()
                gene2 = parts[1].strip()
                weight = float(parts[2])
                edges.append((gene1, gene2, weight))
            except ValueError:
                continue

# ✅ Show top 20 edges
print("Top 20 gene interactions:")
for i, (g1, g2, w) in enumerate(edges[:20]):
    print(f"{i + 1}. {g1} -- {g2} (weight: {w})")

# ✅ Create and draw the graph
G = nx.Graph()
for g1, g2, w in edges:
    G.add_edge(g1, g2, weight=w)

plt.figure(figsize=(10, 8))
pos = nx.spring_layout(G, seed=42)
nx.draw_networkx_nodes(G, pos, node_size=30)
nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.title("Gene Regulatory Network (PIDC)")
plt.axis("off")
plt.tight_layout()
plt.savefig("pidc_graph.png")
plt.show()
