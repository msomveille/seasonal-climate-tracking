#from sys import argv
import csv
import numpy as np
import pandas as pd
import pyemd
from pyemd import emd_with_flow

# Load species seasonal abundance distributions (estimated from eBird data)
abundance_BR = pd.read_table('/Users/marius.somveille/Desktop/niche-tracking/results/seasonalAbundances_B.csv', sep=",")
abundance_NB = pd.read_table('/Users/marius.somveille/Desktop/niche-tracking/results/seasonalAbundances_W.csv', sep=",")

# Load matrix of pairwise distance between every hexagons on the grid
distanceMatrix = np.loadtxt('/Users/marius.somveille/Desktop/niche-tracking/results/distanceMatrix.csv', delimiter=",", skiprows=1)

# Compute optimal redistribution using the Earth Mover's Distance algorithm 
for s in abundance_BR.columns:
  EMD_results = emd_with_flow(np.array(abundance_BR[s]), np.array(abundance_NB[s]), distanceMatrix)
  EMD_results = np.array(EMD_results[1])[(np.array(abundance_BR[s]) > 0),:][:,(np.array(abundance_NB[s]) > 0)]
  np.savetxt("/Users/marius.somveille/Desktop/niche-tracking/results/ORSIM_results_" + s + ".csv", EMD_results, delimiter=',') # Save simulated migratory connectivity

