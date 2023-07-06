# %% Imports

import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms import bipartite
from nxviz import CircosPlot

# %% Import data

dat_raw = pd.read_csv("../../data/Catalog_SiteA.csv", usecols=["LEVEL_ID", "CODE"])

artifact_codes = pd.read_csv("../../data/code_list.csv")

exclude_codes = pd.read_csv("../../data/code_exclude.csv").squeeze("columns")

# Remove bakelite and plastic buttons form exclusions
exclude_codes = exclude_codes[~exclude_codes.isin(["BKLT", "PB"])]

# %% Filter data

dat = dat_raw[~dat_raw["CODE"].isin(exclude_codes)]

# %% Create graph

bpg_assemblage_network = nx.Graph()

bpg_assemblage_network.add_nodes_from(dat["LEVEL_ID"], bipartite="provenience")
bpg_assemblage_network.add_nodes_from(dat["CODE"], bipartite="artifact")

bpg_assemblage_network.add_edges_from(zip(dat["LEVEL_ID"], dat["CODE"]))

# %%

plt.figure()

prov_nodes = nx.bipartite.sets(G=bpg_assemblage_network)[0]

pos = nx.bipartite_layout(bpg_assemblage_network, prov_nodes, align="horizontal")
nx.draw_networkx(
    bpg_assemblage_network, pos=pos, node_size=10, node_color="lightgreen", alpha=0.01
)

plt.show()


# %%
