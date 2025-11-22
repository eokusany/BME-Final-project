#!/usr/bin/env python3
"""
Causal-Learn Example Script
Demonstrates causal discovery using the PC algorithm
"""

import numpy as np
import pandas as pd
from causallearn.search.ConstraintBased.PC import pc
from causallearn.utils.GraphUtils import GraphUtils
import matplotlib.pyplot as plt

# Generate sample data for demonstration
np.random.seed(42)
n_samples = 1000

# Create synthetic causal structure: X -> Y -> Z
X = np.random.normal(0, 1, n_samples)
Y = X + np.random.normal(0, 0.5, n_samples)
Z = Y + np.random.normal(0, 0.5, n_samples)

# Create a DataFrame
data = pd.DataFrame({
    'X': X,
    'Y': Y,
    'Z': Z
})

print("=" * 60)
print("Causal-Learn Example")
print("=" * 60)
print(f"\nSample data shape: {data.shape}")
print(f"\nFirst few rows:")
print(data.head())
print(f"\nData statistics:")
print(data.describe())

# Run PC algorithm for causal discovery
print("\n" + "=" * 60)
print("Running PC Algorithm for Causal Discovery...")
print("=" * 60)

# PC algorithm parameters
cg = pc(data.values, alpha=0.05, indep_test='fisherz', verbose=True)

print("\nCausal graph structure discovered:")
print(cg)

# Visualize the graph
try:
    pdy = GraphUtils.to_pydot(cg)
    pdy.write_png('causal_graph.png')
    print("\n✓ Causal graph saved as 'causal_graph.png'")
except Exception as e:
    print(f"\nNote: Could not save graph visualization: {e}")
    print("(This is okay - the algorithm still ran successfully)")

print("\n" + "=" * 60)
print("Causal discovery completed!")
print("=" * 60)

