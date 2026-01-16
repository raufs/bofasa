# BOFASA Algorithm Documentation

## Table of Contents
1. [Overview](#overview)
2. [Workflow Architecture](#workflow-architecture)
3. [Preparation Phase (`bofasa prep`)](#preparation-phase)
4. [Analysis Phase (`bofasa run`)](#analysis-phase)

---

## Overview

BOFASA (Bacterial Orthology Finding And Syntenic Analysis) implements a multi-stage hierarchical approach to ortholog inference specifically optimized for multi-species bacterial genome analysis. The algorithm prioritizes accuracy over speed, making it ideal for datasets of 4-200 genomes from different species within a genus.

### Key Design Principles

1. **Domain-First Approach**: Proteins are first split into domains, and orthology is determined at the domain level before reconstituting protein-level ortholog groups
2. **Phylogenetic Refinement**: OrthoFinder results are refined using phylogenetic analysis to split paralogs from true orthologs
3. **Syntenic Context Analysis**: Ortholog groups are evaluated based on their genomic neighborhoods
4. **Quality-First Philosophy**: Emphasizes ortholog inference accuracy at the expense of computational speed

---

## Workflow Architecture

BOFASA operates in two distinct phases:

```
┌─────────────────────┐
│   bofasa prep       │  ← Input: Raw genomes (FASTA/GenBank) OR
│   (Preparation)     │           annotation dirs (Bakta/Prokka)
│                     │  → Output: Processed proteins + domains
└──────────┬──────────┘
           │
           ▼
┌─────────────────────┐
│   bofasa run        │  ← Input: Prepared data
│   (Analysis)        │  → Output: Ortholog groups + synteny metrics
└─────────────────────┘
```

---

# Preparation Phase (`bofasa prep`)

The preparation phase processes raw genome data into a format suitable for orthology analysis.

## Step 1. Genome Input Processing

```
Input: FASTA or GenBank files
  ↓
Validation & Format Detection
  ↓
Gene Calling (if FASTA) or CDS Extraction (if GenBank)
```

**Actions:**
- **FASTA files**: Gene calling using Pyrodigal or Prodigal
- **GenBank files**: Direct CDS extraction with optional locus tag renaming
- **Annotation directories**: Can use pre-computed Prokka/Bakta annotations directly
- Locus tags are standardized (default: 3-letter prefix + numeric identifier)

**Key Parameters:**
- `-i, --input-genomes`: Raw genome files in FASTA or GenBank format
- `-a, --annotation-dirs`: Pre-computed Prokka/Bakta annotation directories (alternative to `-i`)
- `-gcm, --gene-calling-method`: Choice between `pyrodigal` (default) or `prodigal`
- `-m, --meta-mode`: Enables meta-genomic gene calling mode (useful for draft genomes or metagenomes)
- `-rlt, --rename-locus-tags`: Forces regeneration of locus tags in GenBank files
- `-rg, --run-genomad`: Run geNomad for annotation of phages and plasmids
- `-emg, --extract-mge-genomes`: Extracts phages/plasmids to use as individual genomes (to use this option, your genomes should be complete)

## Step 2. Protein Extraction
```
Gene Predictions
  ↓
Translate to Proteins
  ↓
Quality Filtering
```

**Actions:**
- Extract protein sequences from CDS features
- Filter incomplete or invalid sequences
- Create sample-specific protein FASTA files

## Step 3. Domain Identification (via PyHMMER+Pfam)
```
Protein Sequences
  ↓
Pfam Domain Annotation
  ↓
Delineation of Proteins by Domain Coordinates
```

**Actions:**
- Use PyHMMER to search for Pfam-A domains in proteins
- Extract domain sequences
- Create "chopped-up" protein FASTA per input genome where each protein is split up into "chunks" based on domain boundary coordinates.

**Reasons:**
- Increases sensitivity by analyzing homology at domain resolution
- Accounts for domain shuffling/fusion/loss events
- Downstream we will be using OrthoFinder to determine course domain ortholog groups, which standardizes for differences in protein-chunk length

## Optional Step 4: Mobile Genetic Element Detection

```
Genomes (if -rg specified)
  ↓
geNomad Analysis for Phage/Plasmid Identification
  ↓
Optional Extraction (if -emg)
```

**Key Parameters:**
- `-rg, --run-genomad`: Enables geNomad for MGE (phage/plasmid) detection
- `-emg, --extract-mge-genomes`: Extracts identified phage/plasmid sequences to be considered as separate genomes

**Effects:**
- `-rg` alone: Annotates MGEs but keeps them in original genomes
- `-rg` + `-emg`: Creates separate genome files for each phage/plasmid

**Rational:**
- Useful to understand how ortholog groups are distributed across autonomous mobile elements
- Can improve resolution of single-copy-core ortholog groups if paralogs exist on MGEs in some genomes

---

# Analysis Phase (`bofasa run`)

The analysis phase performs hierarchical ortholog inference using the prepared domain and protein data.

## Step 1: Inference of Coarse Domain-Resolution Ortholog Groups using OrthoFinder

```
Domain & Inter-domain FASTAs (Bacterial Genomes Only)
  ↓
OrthoFinder (MCL-based clustering)
  ↓
Coarse Domain Ortholog Groups (DOGs)
```

**Key Parameters for OrthoFinder:**
- `-mi, --mcl-inflation` (default: 1.2): MCL inflation parameter
  - Lower values (0.8-1.2): More coarse clustering, larger ortholog groups
  - Higher values (1.2-5.0): Finer clustering, more granular groups
- `-us, --ultra-sens`: Use DIAMOND ultra-sensitive mode
  - Recommended for highly divergent species

> [!NOTE]
> When MGEs are present (from `bofasa prep -emg`), they are **excluded** from OrthoFinder to prevent biasing the core ortholog group structure. MGEs are integrated in Step 1b.

## Step 2: MGE Integration (_auxiliary_; if `-emg` used in prep)

```
MGE Protein Chunks
  ↓
Reflexive DIAMOND blastp Alignment (MGE protein chunks vs. MGE protein chunks)
  ↓
Single-Linkage Clustering (using slclust by Brian Haas)
  ↓
MGE Homology Clusters
  ↓
DIAMOND Alignment (MGE protein chunks vs bacterial chromosome protein chunks)
  ↓
OG Assignment (best hit or singleton)
  ↓
Construction of Modified  Domain-Resolution Coarse Ortholog Group by Genome/MGE Matrix
```

**Algorithm:**
1. **Concatenate MGE protein chunks**: All MGE protein chunks into single FASTA
2. **Concatenate bacterial protein chunks**: All bacterial chromosome protein chunks into single FASTA
3. **Create DIAMOND databases**: Separate databases for MGE and bacterial proteins
4. **Reflexive MGE alignment**: DIAMOND blastp of MGEs against themselves
   - Filter hits: require ≥50% query coverage AND ≥50% subject coverage AND E-value < 1e-3
   - Purpose: Identify MGE protein chunks that are homologous to each other
5. **Cluster MGE protein chunks**: Single-linkage clustering (using slclust) on filtered pairs are used to groups them into homologous clusters
6. **MGE vs chromosomal protein chunk alignment**: DIAMOND blastp of MGE protein chunks against bacterial chromosome protein chunks database
   - Sorted by bitscore (descending) using the commandline program `sort`
   - Each MGE protein chunk assigned to OG of its best chromosomal hit
7. **Cluster consensus**: For each MGE homology cluster:
   - If all members map to same OG → assign all to that OG
   - If members map to different OGs → leave OG assignments based on best hits for individual protein chunks and leave the rest of the protein chunks in the MGE homologous cluster as singletons
8. **Output modified tables**: Modified files can be found under the `OrthoFinder_Results/` directory.
   - `Orthogroups_Modified.tsv`: Original OGs + MGE columns (includes some previously single (aka "unassigned") protein chunks that are homologous to protein chunks from MGEs)
   - `Orthogroups_UnassignedGenes_Modified.tsv`: Singletons + new MGE singletons

## Step 3: Phylogenetic Refinement of Domain-Resolution Ortholog Groups

```
Coarse Domain-Resolution Ortholog Groups from OrthoFinder
  ↓
Per-Group Multiple Sequence Alignment using MUSCLE super5
  ↓
Phylogenetic Tree Construction using FastTree2
  ↓
Tree Analysis to Split Paralogs
  ↓
Refined Domain Ortholog Groups
```

**Algorithm:**
For each **coarse domain-resolution ortholog group**:
   - Build multiple sequence alignment (default: MUSCLE).
   - Construct phylogenetic tree (default: FastTree2).
   - Apply midpoint-rooting. If `-rs N` flag is specified, then root the tree randomly _N_-1 amount of times and once using midpoint. 
   - Partition phylogeny into **refined domain-resolution ortholog groups** using using bofasa's recursive algorithm.
   - Further refine/adjust delineations of refined ortholog groups: (1) further split up disjoint ortholog groups [*see note* :scissors:] and (2) fixation index calculation to assess whether split partitions should be merged back [*see note* :tent:].

### Overview of bofasa's recursive algorithm:

Goal: Split a phylogenetic tree into ortholog groups by minimizing gene duplications per sample

#### Scoring System:

- ***Tree score*** = $\sum_{i=1}^{n} |(x_i)-1|$, where _n_ is the number of genomes/samples and $x_i$ is the copy count of proteins from genome _i_ at the particular node being assessed of the ortholog group phylogeny/tree.
   - Score of 0 = perfect single-copy ortholog group (one gene per sample)
   - Higher scores = more duplications/paralogs present 
   - Also considers total branch length as a secondary criterion

#### Return: List of partitions (sets of sequences)

#### Recursive Logic:
- Base cases (stop splitting):
   - For proteins: if total branch length = 0
   - For domains: if branch length ≤ minimum threshold
   - If only one sample present
   - If no split reduces the score

- To find the best split(s) each iteration:
   - Calculate score for entire tree
   - Calculate score for every subtree (traverse preorder)
   - Identify subtree(s) with minimum score
   - Select non-overlapping subtrees that minimize duplication

#### Handle remaining sequences:
- Create "complement" partition (sequences not in best subtrees)
- ***If complement has ≥2 samples, recursively split it***
- Combine best subtrees + complement splits
  
### Details on addition refinement steps:

#### :scissors: Further split up disjoint ortholog groups

Ensures each final ortholog group form clear evolutionary units (monophyletic):

```
For each split group G:
  if G forms a monophyletic clade:
    accept G
  else:
    # G is paraphyletic or polyphyletic
    # Split into largest monophyletic sub-clades
    extract maximal monophyletic subsets
```

This prevents ortholog groups from spanning multiple independent evolutionary lineages.

#### :tent: Fixation index calculation of innernodes of coarse ortholog group phylogeny to assess whether split partitions should be merged back 

Re-assesses each internal node of coarse ortholog group phylogenies that include multiple refined ortholog groups to see if they are evolutionarily more appropriate to group together.

- Uses largest-to-smallest node traversion
- Calculates fixation index (FST) between split groups
- If FST < threshold → merge groups back together (default threshold = 0.1)
- Prevents over-splitting of true orthologs

**Key Parameters:**
- `-spr, --skip-phylo-refine`: Skip phylogenetic refinement entirely
  - Uses OrthoFinder results directly
  - Much faster but less accurate
  
- `-fic`, `--fixation-index-cutoff` (default: 0.1): Minimum FST to accept a split
  - Lower values (0.1-0.24): More re-merging
  - Higher values (0.26-1.0): Less re-merging
  
- `-smb, --skip-merge-back`: Skip merge-back assessment
  - Faster but may over-split ortholog groups
  
- `-rs, --rooting-seeds` (default: 1): Number of outgroup rooting attempts
  - Higher values: More robust rooting but slower
  - 1 (default): Fastest, uses midpoint rooting
  - 100: Tries midpoint rooting + 99 random nodes as roots to see if they result in higher quality refined ortholog groups
  - Best partioning is selected based on scoring:
     - Sort by: (1) tree score, (2) number of groups, (3) total branch length and choose partition with lowest combined score
  
- `-qa, --quality-alignments`: Use high-quality alignment settings
  - MAFFT instead of MUSCLE
  - IQ-TREE instead of FastTree2
  - Much slower but should lead to higher quality results

## Step 4: Protein Ortholog Group Determination

```
Refined Domain Ortholog Groups
  ↓
Protein-to-Domain Mapping
  ↓
Jaccard Similarity Calculation Between Proteins Based on Domain Ortholog Groups Composition
  ↓
MCL (Graph-Based) Protein Clustering
  ↓
Coarse Protein Ortholog Groups
  ↓
Distance Based Tree Construction using FastME
  ↓
Tree Refinement (for multi-copy groups)
  ↓
Final Protein Ortholog Groups
```

### Part A: Initial Jaccard-Based Clustering
1. For each protein, identify its domain ortholog group memberships
2. Calculate Jaccard similarity between proteins:
   ```
   J(A,B) = |DOGs(A) ∩ DOGs(B)| / |DOGs(A) ∪ DOGs(B)|
   ```
3. Create protein similarity graph (edges = similarity ≥ threshold; default threshold is )
4. Identify connected components = coarse protein ortholog groups

### Part B: Phylogenetic Refinement of Multi-Copy Groups
For protein ortholog groups meeting refinement criteria (≥2 copies in any genome, ≥2 genomes, ≥4 proteins total):

1. **Domain-based distance calculation**:
   - Represent each protein as a vector of domain ortholog group counts
   - Calculate pairwise cosine distances between protein domain vectors
   - Create distance matrix in PHYLIP format

2. **Coarse protein ortholog group tree construction**:
   - Build phylogenetic tree using FastME

3. **Tree-based splitting**:
   - Root tree at midpoint
   - Apply recursive splitting at duplication nodes (uses same recusive function as for domain-resolution OG splitting)
       - Note, however, that there is a slight difference for domain-resolution and protein-resolution recursive splitting. Namely, for proteins, we stop only if the branch length sum = 0 (identical) whereas for domains we stop if branch length ≤ threshold. This is inorder to account for pseudovalues added by FastTree 2.          

4. **Outlier artifact removal**:
   - For each protein in a split group:
     - Calculate max Jaccard similarity to proteins within the group
     - Calculate max Jaccard similarity to proteins outside the group
     - If external similarity ≥ internal similarity → remove as singleton
   - Prevents misplaced proteins from remaining in inappropriate groups - could result from nesting pattern in FastME tree that arises akin to similar artifacts in neighbor-joining trees for instance.

> [!NOTE]
> When calculating Jaccard similarity between proteins based on domain-resolution ortholog groups (DOGs) for outlier detection, DOGs that are largely (>80%) inter-domain regions (IDR) are included to provide higher resolution. However, for initial aggregation of proteins into coarse protein ortholog groups, the Jaccard similarity indices measured does not account for DOGs as they might contribute noise.

5. **Monophyletic enforcement**:
   - Check if each split group forms a monophyletic clade
   - If non-monophyletic, further split into largest monophyletic sub-clades
   - Ensures evolutionary coherence within final ortholog groups

**Key Parameters:**
- `-dj, --dog-jaccard` (default: 0.5): Jaccard similarity threshold for initial clustering
  - Lower values (0.1-0.49): More permissive, allows for more domain composition differences
  - Higher values (0.51-1.0): Stricter, requires domain compositions between proteins to be high

## Step 5: Syntenic Context Analysis

```
Protein Ortholog Groups
  ↓
Extract Genomic Neighborhoods
  ↓
Calculate Context Conservation Score
  ↓
Ortholog Groups + Synteny Metrics
```

### Context Conservation Score

#### Formula:

$Context\ Conservation\ Score = Avg.\ Neighbor\ Count/Total\ Distinct\ Neighbors$

- **Avg. Neighbor Count**: The average of the number of unique OGs that appear as neighbors to the focal OG ***per*** context.
- **Total Distinct Neighbors**: Total number of unique OGs that ever appear as neighbors to the focal OG across ***gall*** context.

#### Interpretation:
- Score → 1: High conservation (same neighbors always present)
- Score → 0: Low conservation (many different neighbors, the focal ortholog group is found in very different contexts)

#### Example:
- OG appears 10 times
- Has 5 unique neighbor OGs across all instances
- Average of 3 neighbors per instance
Score = 3/5 = 0.6 (moderately conserved context)

#### Two versions:
- Standard: All contexts (including near scaffold/contig edges)
- Complete: Only complete contexts (away from scaffold/contig edges)

#### Key Parameters:
- `-sr, --surrounding-bp` (default: 10000): Base pairs to analyze around each gene

## Step 6: Final Report Generation

```
All Analysis Results
  ↓
Data Integration
  ↓
Excel Spreadsheet Generation & HTML Visualization
```

### Visualization Overview:
- Interactive HTML plot showing:
  - X-axis: Ortholog group conservation (conservation across genomes)
  - Y-axis: Context conservation score (syntenic conservation)
  - ***Hover***: Detailed OG information

## Optional Step 7: Consensus Sequence + Profile HMM Generation for Protein Ortholog Groups

### Premise:

You can also use the option `-ogc` option to generate a multi-FASTA containing consensus sequences of each ortholog groups as well as profile HMMs for them. 

### Key Parameter:
- `-ogc, --og-consensus`: Generate consensus sequences for each ortholog group

### Algorithm:
1. Align all sequences in each ortholog group
2. Determine consensus sequence + profile HMM for each ortholog group

## Optional Step 8: Core Genome Alignment

### Premise:

Construct a concatenated multi-FASTA alignment of the strict or loose single-copy-core genome for downstream phylogenomics.

### Key Parameter:
- `-cg, --core-genome`: Construct concatenated core genome alignment
- `-ns, --near-scc-prop` (default: 0.95): Minimum proportion of genomes for an ortholog group to be considered part of the single-copy-core.
  - Lower values (0.80-0.90): More permissive, more genes
  - Higher values (0.98-1.00): Stricter core, fewer genes

> [!IMPORTANT]
> An ortholog group can still be considered part of the loose single-copy-core (`-ns` <1.0) when genomes have multiple copies of it. These are ignored similar to the absence of the gene in other genomes. For instance, if your dataset has 100 genomes and 97 of the genomes have the focal ortholog group in single-copy, but 2 genomes lack it and 1 genome has two copies of the ortholog group, if `-ns` >= 0.97, then it is still considered part of the scc and treated as absent in that 1 genome with paralogs. **Also, if you are _not_ using a strict core genome (`-ns` set to 1.0), we recommend that you use a partition-based approach for phylogeny modeling and using the individual ortholog group alignments folder as input to IQ-TREE instead of the concatenated multi-FASTA file..**
  
### Algorithm:
1. Identify single-copy core ortholog groups (SCC-OGs)
   - Present in ≥95% of genomes by default unless `-ns` modified.
2. Align each SCC-OG
3. Concatenate alignments and output multi-FASTA suitable for phylogenomics

## Key Citations

When using BOFASA, please cite:

> Salamzade, R., Kottapalli, A., & Kalan, L. (2025). bofasa: high-quality orthology inference across multiple bacterial species.

And the underlying tools:
- **OrthoFinder**: Emms, D.M. & Kelly, S. (2019) Genome Biology
- **geNomad**: Camargo, A.P. et al. (2023) Nature Biotechnology
- **FastTree**: Price, M.N. et al. (2010) PLoS One (domain-level phylogeny)
- **FastME**: Lefort, V. et al. (2015) Molecular Biology and Evolution (protein-level phylogeny)
- **MUSCLE**: Edgar, R.C. (2022) Nature Communications
