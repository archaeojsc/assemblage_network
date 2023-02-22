# %% Imports

import pandas as pd
import networkx as nx

# %% Import data

dat_raw = pd.read_csv("..\data\Catalog_SiteA.csv", usecols=["LEVEL_ID", "CODE"])

artifact_codes = pd.read_csv("..\data\code_list.csv")

exclude_codes = pd.read_csv("../data/code_exclude.csv").squeeze("columns")


# %% Filter data
