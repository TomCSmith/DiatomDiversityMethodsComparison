import pandas as pd
import numpy as np
import os

"""
Description
-----------
Calculating the Species Richness, Total Number of Species, Species Diversity,
and SAD fits (dgamma, lognormal, poisson lognormal, ZSM (?), logseries).

"""

#Specify datapaths
holmes_data_path = os.path.join("..", "data", "formatted",
        "HolmesDiatomsAbundByLake_formatted.csv")

lake_data_path = os.path.join("..", "data", "formatted",
    "LakeData_formatted.csv")

comm_data = pd.read_csv(holmes_data_path)
lake_data = pd.read_csv(lake_data_path)

comm_data.set_index("LakeName", inplace=True)
lake_data.set_index("LakeName", inplace=True)

# Total community abundance
N_matrix = pd.DataFrame(comm_data.sum(axis=1), columns=["N"])

# Total species richness per community
ind_matrix = comm_data != 0
S_matrix = pd.DataFrame(ind_matrix.sum(axis=1), columns=["S"])

# Shannon Species diversity
# TODO

# FIT SADs
# TODO

N_and_S_data = lake_data.join([N_matrix, S_matrix])

N_and_S_data.to_csv(os.path.join("..", "data", "formatted",
    "Lake_SAD_analysis.csv"))
