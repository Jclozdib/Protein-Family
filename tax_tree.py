import networkx as nx
import matplotlib.pyplot as plt
from collections import defaultdict
from networkx.drawing.nx_agraph import graphviz_layout
import random 

input_file = "lineage.txt"

with open(input_file, "r") as f:
    data = [line.strip() for line in f]

lineages = []
for entry in data:
    _, lineage_str = entry.split("\t")
    lineages.append(lineage_str.split(" -> "))

abundance = defaultdict(int)
for lineage in lineages:
    for group in lineage:
        abundance[group] += 1

tree = nx.DiGraph()
for lineage in lineages:
    for i in range(len(lineage) - 1):
        parent = lineage[i]
        child = lineage[i + 1]
        tree.add_edge(parent, child)

node_sizes = [abundance[node] * 100 for node in tree.nodes]  # Scaling 
color_list = ['red', 'blue', 'green', 'orange', 'purple', 'yellow', 'cyan', 'magenta']

# Assigning random color just for a better visualization
color_map = {node: random.choice(color_list) for node in tree.nodes}
node_colors = [color_map[node] for node in tree.nodes]

pos = graphviz_layout(tree, prog="dot")

y_offset = -0.1
x_offset = 0.1

label_pos = {node: (pos[node][0] + x_offset, pos[node][1] + y_offset) for node in tree.nodes}

node_sizes = [abundance[node] * 100 for node in tree.nodes]  # Scale node sizes

plt.figure(figsize=(16, 12))

nx.draw(
    tree, pos, with_labels=False, 
    node_size=node_sizes, node_color=node_colors, 
    font_size=8, font_weight="bold"
)

for node, (x, y) in label_pos.items():
    plt.text(
        x, y, node, fontsize=8, fontweight='bold', color='black', 
        ha='center', va='center', rotation=20  
    )

plt.title("Taxonomic Tree with Node Sizes Proportional to Relative Abundance", fontsize=16)
plt.axis("off")  

output_file = "taxonomic_tree.png"
plt.savefig(output_file) 

print(f"Taxonomic tree saved to {output_file}")
