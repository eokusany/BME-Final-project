#!/usr/bin/env python3
"""
Causal Graph Discovery using PC Algorithm
Implements constraint-based causal discovery to infer causal structure

Methods:
- PC Algorithm: Constraint-based causal discovery using causal-learn library
- MR-based graph: Uses MR results to create directed gene -> outcome edges
- Bootstrap sampling: Creates multiple observations for PC algorithm

Development notes:
- Initially tried using SNP-level data directly but PC algorithm found 0 edges
- Switched to bootstrap sampling approach after reading PC algorithm documentation
- MR-based graph provides complementary view using actual MR estimates

Author: Emmanuel Okusanya
Date: 2024
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from sklearn.preprocessing import StandardScaler
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Import causal-learn for PC algorithm
# Using causal-learn library - found this after trying to implement PC from scratch
# Much easier to use existing implementation!
try:
    from causallearn.search.ConstraintBased.PC import pc
    from causallearn.utils.GraphUtils import GraphUtils
    CAUSAL_LEARN_AVAILABLE = True
except ImportError:
    CAUSAL_LEARN_AVAILABLE = False
    print("Warning: causal-learn not available. Please install: pip install causal-learn")

print("=" * 70)
print("Improved Causal Graph Discovery")
print("=" * 70)
print()

# ========= 1. LOAD DATA =========

print("Loading MR data...")
exposure = pd.read_csv("exposure.txt", sep="\t")
outcome = pd.read_csv("outcome.txt", sep="\t")

print(f"  - Exposure SNPs: {len(exposure):,}")
print(f"  - Outcome SNPs: {len(outcome):,}")

# MERGE AND HARMONIZE (same as run.py) 

print("\nMerging and harmonizing data...")

# Rename columns into a standard format 
exposure = exposure.rename(columns={
    "rsid": "snp",
    "effect_allele": "ea_exposure",
    "other_allele": "oa_exposure",
    "beta": "beta_exposure",
    "se": "se_exposure",
    "pval": "pval_exposure",
    "n": "n_exposure"
})

# Preserve gene_name if present
has_gene_name = "gene_name" in exposure.columns
if has_gene_name:
    exposure["gene_name"] = exposure["gene_name"]

outcome = outcome.rename(columns={
    "rsid": "snp",
    "effect_allele": "ea_outcome",
    "other_allele": "oa_outcome",
    "beta": "beta_outcome",
    "se": "se_outcome",
    "pval": "pval_outcome",
    "n": "n_outcome"
})

# Merge on SNP
merged = exposure.merge(outcome, on="snp", how="inner")
print(f"  - Merged SNPs: {merged.shape[0]:,}")

# Define palindromic pairs
palindromic_pairs = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}

def is_palindromic(a1, a2):
    return (a1, a2) in palindromic_pairs

def harmonise_row(row):
    ea_exp, oa_exp = row["ea_exposure"], row["oa_exposure"]
    ea_out, oa_out = row["ea_outcome"], row["oa_outcome"]
    beta_out = row["beta_outcome"]

    # Drop palindromic SNPs
    if is_palindromic(ea_exp, oa_exp):
        return pd.Series({"keep": False})

    # Alleles already aligned
    if (ea_exp == ea_out) and (oa_exp == oa_out):
        return pd.Series({
            "beta_out_harmonised": beta_out,
            "keep": True
        })

    # Alleles reversed -> flip outcome beta
    if (ea_exp == oa_out) and (oa_exp == ea_out):
        return pd.Series({
            "beta_out_harmonised": -beta_out,
            "keep": True
        })

    # Anything else -> drop
    return pd.Series({"keep": False})

harm = merged.apply(harmonise_row, axis=1)
merged_h = pd.concat([merged, harm], axis=1)
merged_h = merged_h[merged_h["keep"]].copy()
merged_h["beta_outcome"] = merged_h["beta_out_harmonised"]

print(f"  - SNPs after harmonization: {merged_h.shape[0]:,}")

# Filter out very small beta_exposure
min_beta_exp = 1e-6
merged_h = merged_h[np.abs(merged_h["beta_exposure"]) > min_beta_exp].copy()
print(f"  - SNPs after beta filter: {merged_h.shape[0]:,}")

# ========= 3. CALCULATE GENE-LEVEL SUMMARY STATISTICS =========

print("\nCalculating gene-level summary statistics...")

genes = sorted(merged_h['gene_name'].unique())
print(f"  - Genes: {len(genes)}")

# Calculate gene-level IVW estimates (same as run.py)
gene_summary = []

for gene in genes:
    gene_data = merged_h[merged_h['gene_name'] == gene].copy()
    
    if len(gene_data) == 0:
        continue
    
    # Calculate Wald ratios
    beta_x = gene_data["beta_exposure"].values
    beta_y = gene_data["beta_outcome"].values
    se_x = gene_data["se_exposure"].values
    se_y = gene_data["se_outcome"].values
    
    # Wald ratio for each SNP
    beta_mr = beta_y / beta_x
    
    # Delta-method SE
    se_mr = np.sqrt(
        (se_y ** 2) / (beta_x ** 2) +
        (beta_y ** 2) * (se_x ** 2) / (beta_x ** 4)
    )
    
    # IVW MR for this gene
    w = 1.0 / (se_mr ** 2)
    beta_ivw = np.sum(w * beta_mr) / np.sum(w)
    se_ivw = np.sqrt(1.0 / np.sum(w))
    z_ivw = beta_ivw / se_ivw
    
    # P-value calculation
    if np.abs(z_ivw) > 6:
        log_p = stats.norm.logsf(np.abs(z_ivw)) + np.log(2)
        p_ivw = np.exp(log_p)
    else:
        p_ivw = 2 * (1 - stats.norm.cdf(np.abs(z_ivw)))
    
    # Gene-level exposure and outcome effects
    # Use IVW-weighted average for exposure
    w_exp = 1.0 / (se_x ** 2)
    beta_exp_gene = np.sum(w_exp * beta_x) / np.sum(w_exp)
    se_exp_gene = np.sqrt(1.0 / np.sum(w_exp))
    
    # Use IVW-weighted average for outcome
    w_out = 1.0 / (se_y ** 2)
    beta_out_gene = np.sum(w_out * beta_y) / np.sum(w_out)
    se_out_gene = np.sqrt(1.0 / np.sum(w_out))
    
    gene_summary.append({
        'gene': gene,
        'n_snps': len(gene_data),
        'beta_exp': beta_exp_gene,
        'se_exp': se_exp_gene,
        'beta_out': beta_out_gene,
        'se_out': se_out_gene,
        'beta_mr_ivw': beta_ivw,
        'se_mr_ivw': se_ivw,
        'pval_mr': p_ivw,
        'z_mr': z_ivw
    })

gene_df = pd.DataFrame(gene_summary)
print(f"  - Gene-level summary calculated for {len(gene_df)} genes")

# ========= 4. PREPARE DATA FOR PC ALGORITHM =========

print("\nPreparing data for PC algorithm...")


n_bootstrap = 1000  # Number of bootstrap samples
# Chose 1000 as a balance - enough for stable results but not too slow
# Initially tried 100 but PC algorithm found very few edges
# Increased to 1000 after reading that PC algorithm needs sufficient sample size
print(f"  - Creating {n_bootstrap} bootstrap samples for PC algorithm...")

np.random.seed(42)  # For reproducibility - learned this is important!
bootstrap_data = []

for _ in range(n_bootstrap):
    sample = {}
    
    # Sample SNPs for each gene
    for gene in genes:
        gene_snps = merged_h[merged_h['gene_name'] == gene]
        if len(gene_snps) > 0:
            # Randomly sample a SNP from this gene
            sampled_snp = gene_snps.sample(n=1, replace=True).iloc[0]
            sample[gene] = sampled_snp['beta_exposure']
        else:
            sample[gene] = 0.0
    
    # Sample outcome (use average across all SNPs)
    sample['LVEF'] = merged_h.sample(n=min(100, len(merged_h)), replace=True)['beta_outcome'].mean()
    
    bootstrap_data.append(sample)

data_df = pd.DataFrame(bootstrap_data)
print(f"  - Data matrix shape: {data_df.shape}")

# Standardize the data
scaler = StandardScaler()
data_scaled = scaler.fit_transform(data_df.values)

var_names = list(genes) + ['LVEF']

# ========= 5. RUN PC ALGORITHM =========
# PC algorithm: Constraint-based causal discovery
# Tests conditional independence to infer causal structure
# Alpha = significance threshold for independence tests
# Using Fisher's Z test (appropriate for continuous data)

G_pc = None
if CAUSAL_LEARN_AVAILABLE:
    print("\nRunning PC algorithm...")
    
    alpha = 0.05  # Standard significance level - could adjust but 0.05 is conventional
    
    try:
        cg = pc(data_scaled, alpha=alpha, indep_test='fisherz', stable=True)
        
        print(f"  - PC algorithm completed")
        print(f"  - Number of Edges: {cg.G.get_num_edges()}")
        print(f"  - Number of Nodes: {cg.G.get_num_nodes()}")
        
        # Create NetworkX graph from PC result
        G_pc = nx.Graph()
        G_pc.add_nodes_from(var_names)
        
        # Add edges from PC graph
        for i in range(cg.G.get_num_nodes()):
            for j in range(i+1, cg.G.get_num_nodes()):
                if cg.G.graph[i, j] != 0 or cg.G.graph[j, i] != 0:
                    G_pc.add_edge(var_names[i], var_names[j])
        
        print(f"  - PC graph edges: {G_pc.number_of_edges()}")
        
    except Exception as e:
        print(f"  - Error running PC algorithm: {e}")
        print("  - Will create MR-based graph...")

# ========= 6. CREATE MR-BASED CAUSAL GRAPH =========

print("\nCreating MR-based causal graph...")

G_mr = nx.DiGraph()
G_mr.add_nodes_from(var_names)

# Add edges based on MR results
# Gene -> LVEF edges based on MR significance
pval_threshold = 0.05

for _, row in gene_df.iterrows():
    gene = row['gene']
    pval = row['pval_mr']
    beta_mr = row['beta_mr_ivw']
    
    if pval < pval_threshold:
        # Add directed edge from gene to LVEF
        G_mr.add_edge(gene, 'LVEF',
                     weight=abs(beta_mr),
                     beta_mr=beta_mr,
                     pval=pval,
                     n_snps=row['n_snps'])

print(f"  - MR graph nodes: {G_mr.number_of_nodes()}")
print(f"  - MR graph edges: {G_mr.number_of_edges()}")

# Add gene-gene edges based on correlation
# Calculate correlation between gene expression effects
gene_corr_data = []
for gene in genes:
    gene_snps = merged_h[merged_h['gene_name'] == gene]
    if len(gene_snps) > 0:
        # Use average beta_exp as gene expression level
        gene_corr_data.append({
            'gene': gene,
            'beta_exp': gene_snps['beta_exposure'].mean()
        })

if len(gene_corr_data) > 1:
    gene_corr_df = pd.DataFrame(gene_corr_data)
    # Calculate pairwise correlations using SNP-level data
    gene_pairs = []
    for i, gene1 in enumerate(genes):
        for gene2 in genes[i+1:]:
            snps1 = merged_h[merged_h['gene_name'] == gene1]
            snps2 = merged_h[merged_h['gene_name'] == gene2]
            
            # Find common SNPs or use all SNPs
            # For simplicity, use correlation of average effects
            if len(snps1) > 0 and len(snps2) > 0:
                # Sample-based correlation
                n_samples = min(100, len(snps1), len(snps2))
                sample1 = snps1.sample(n=n_samples, replace=True)['beta_exposure'].values
                sample2 = snps2.sample(n=n_samples, replace=True)['beta_exposure'].values
                corr = np.corrcoef(sample1, sample2)[0, 1]
                
                if not np.isnan(corr) and abs(corr) > 0.1:  # Threshold
                    G_mr.add_edge(gene1, gene2, weight=abs(corr), correlation=corr)

# ========= 7. VISUALIZE GRAPHS =========

print("\nVisualizing causal graphs...")

# Create radial layout
def create_radial_layout(genes, center_node='LVEF'):
    pos = {}
    pos[center_node] = (0, 0)
    
    n_genes = len(genes)
    radius = 3.0
    for i, gene in enumerate(genes):
        angle = 2 * np.pi * i / n_genes
        x = radius * np.cos(angle)
        y = radius * np.sin(angle)
        pos[gene] = (x, y)
    
    return pos

# Visualize MR-based graph
plt.figure(figsize=(18, 18))
pos = create_radial_layout(genes, 'LVEF')

# Draw edges
edges_mr = G_mr.edges()
if len(edges_mr) > 0:
    weights_mr = [G_mr[u][v].get('weight', 1.0) for u, v in edges_mr]
    max_weight = max(weights_mr) if weights_mr else 1.0
    edge_widths = [w * 3 / max_weight for w in weights_mr]
    
    nx.draw_networkx_edges(G_mr, pos, edgelist=edges_mr,
                          width=edge_widths,
                          alpha=0.6,
                          edge_color='gray',
                          arrows=True,
                          arrowsize=25,
                          arrowstyle='->')

# Draw nodes
gene_nodes = [n for n in G_mr.nodes() if n != 'LVEF']
nx.draw_networkx_nodes(G_mr, pos, nodelist=gene_nodes,
                      node_color='lightblue',
                      node_size=2500,
                      alpha=0.8,
                      edgecolors='black',
                      linewidths=2)

nx.draw_networkx_nodes(G_mr, pos, nodelist=['LVEF'],
                      node_color='orange',
                      node_size=3500,
                      alpha=0.9,
                      edgecolors='black',
                      linewidths=3)

# Labels
labels = {n: n for n in G_mr.nodes()}
nx.draw_networkx_labels(G_mr, pos, labels,
                       font_size=10,
                       font_weight='bold',
                       font_color='black')

plt.title("MR-Based Causal Graph\n(Gene Expression → LVEF, weighted by MR effect size)",
          fontsize=18, fontweight='bold', pad=20)
plt.axis('off')
plt.tight_layout()
plt.savefig('causal_graph_mr_improved.png', dpi=300, bbox_inches='tight', facecolor='white')
print("  - Saved MR graph to: causal_graph_mr_improved.png")

# Visualize PC algorithm graph if available
if G_pc is not None and G_pc.number_of_edges() > 0:
    plt.figure(figsize=(18, 18))
    pos_pc = create_radial_layout(genes, 'LVEF')
    
    nx.draw_networkx_edges(G_pc, pos_pc,
                          alpha=0.5,
                          edge_color='blue',
                          width=2,
                          style='solid')
    
    gene_nodes_pc = [n for n in G_pc.nodes() if n != 'LVEF']
    nx.draw_networkx_nodes(G_pc, pos_pc, nodelist=gene_nodes_pc,
                          node_color='lightgreen',
                          node_size=2500,
                          alpha=0.8,
                          edgecolors='black',
                          linewidths=2)
    
    nx.draw_networkx_nodes(G_pc, pos_pc, nodelist=['LVEF'],
                          node_color='orange',
                          node_size=3500,
                          alpha=0.9,
                          edgecolors='black',
                          linewidths=3)
    
    labels_pc = {n: n for n in G_pc.nodes()}
    nx.draw_networkx_labels(G_pc, pos_pc, labels_pc,
                           font_size=10,
                           font_weight='bold',
                           font_color='black')
    
    plt.title("PC Algorithm: Discovered Causal Structure\n(Constraint-based causal discovery)",
              fontsize=18, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig('causal_graph_pc_improved.png', dpi=300, bbox_inches='tight', facecolor='white')
    print("  - Saved PC graph to: causal_graph_pc_improved.png")

# ========= 8. SAVE GRAPH DATA =========

print("\nSaving graph data...")

# Save MR graph edges
mr_edges_list = []
for u, v, data in G_mr.edges(data=True):
    mr_edges_list.append({
        'source': u,
        'target': v,
        'weight': data.get('weight', 1.0),
        'beta_mr': data.get('beta_mr', np.nan),
        'pval': data.get('pval', np.nan),
        'correlation': data.get('correlation', np.nan),
        'n_snps': data.get('n_snps', np.nan)
    })

pd.DataFrame(mr_edges_list).to_csv('causal_graph_mr_edges_improved.csv', index=False)
print("  - Saved MR edges to: causal_graph_mr_edges_improved.csv")

# Save gene summary statistics
gene_df.to_csv('gene_mr_summary.csv', index=False)
print("  - Saved gene MR summary to: gene_mr_summary.csv")

# Save PC graph edges if available
if G_pc is not None:
    pc_edges_list = []
    for u, v in G_pc.edges():
        pc_edges_list.append({'source': u, 'target': v})
    pd.DataFrame(pc_edges_list).to_csv('causal_graph_pc_edges_improved.csv', index=False)
    print("  - Saved PC edges to: causal_graph_pc_edges_improved.csv")


print("\n" + "=" * 70)
print("Summary")
print("=" * 70)
print(f"Total genes: {len(genes)}")
print(f"MR graph: {G_mr.number_of_nodes()} nodes, {G_mr.number_of_edges()} edges")
if G_pc is not None:
    print(f"PC graph: {G_pc.number_of_nodes()} nodes, {G_pc.number_of_edges()} edges")

print(f"\nGenes significantly associated with LVEF (p < 0.05):")
sig_genes = gene_df[gene_df['pval_mr'] < 0.05].sort_values('pval_mr')
for _, row in sig_genes.iterrows():
    print(f"  {row['gene']}: β = {row['beta_mr_ivw']:.4f}, p = {row['pval_mr']:.3e}, n_SNPs = {row['n_snps']}")

print(f"\nGenes connected to LVEF in MR graph: {G_mr.degree('LVEF')}")

print("\n" + "=" * 70)
print("Processing complete!")
print("=" * 70)

