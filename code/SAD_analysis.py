from __future__ import division
import pandas as pd
import numpy as np
import os

# Personal Package
import macroeco.models as mod
import macroeco.empirical as emp

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
SD_vec = []
for index in comm_data.index:

    sad = comm_data.ix[index][comm_data.ix[index] != 0]
    N = np.sum(sad)
    tSD = -np.sum((sad / N) * np.log(sad / N))
    SD_vec.append(tSD)

SD_matrix = pd.DataFrame(SD_vec, columns=["Shannon_Diversity"],
    index=S_matrix.index)

# FIT SADs
lognorm_params = []
logser_params = []
dgamma_params = []
for index in comm_data.index:

    print("Fitting SADs...")
    sad = comm_data.ix[index][comm_data.ix[index] != 0]
    lognorm_params.append(mod.lognorm.fit_mle(sad, fix_mean=True))
    logser_params.append(mod.logser_uptrunc.fit_mle(sad))
    dgamma_params.append(mod.dgamma.fit_mle(sad))

# Unpack SADs
mus, sigmas = zip(*lognorm_params)
ps, bs = zip(*logser_params)
alphas, thetas, bs = zip(*dgamma_params)
fishers_pseudo_alpha = S_matrix.S * (-1 / np.log(1 - np.array(ps)))

params_matrix = pd.DataFrame(np.array([sigmas, ps, fishers_pseudo_alpha,
        alphas, thetas]).T, columns=["sigma", "p", "fishers_alpha",
        "alpha", "theta"], index=S_matrix.index)

# Make and output analysis
N_and_S_data = lake_data.join([N_matrix, S_matrix, SD_matrix, params_matrix])

# Drop Dana because it doesn't have a community vector
N_and_S_data = N_and_S_data.drop("Dana")

N_and_S_data.to_csv(os.path.join("..", "data", "formatted",
    "Lake_SAD_analysis.csv"))
