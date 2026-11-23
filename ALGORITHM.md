# BOFASA Algorithm Documentation

## Table of Contents
1. [Overview](#overview)
2. [Workflow Architecture](#workflow-architecture)
3. [Preparation Phase (`bofasa prep`)](#preparation-phase)
4. [Analysis Phase (`bofasa run`)](#analysis-phase)
5. [Database Setup (`bofasa setup`)](#database-setup-bofasa-setup)
6. [Standalone Pfam Domain Matrix (`bofasa domscape`)](#standalone-pfam-domain-matrix-bofasa-domscape)
7. [Reproducibility](#reproducibility)

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

After a one-time `bofasa setup` to download databases, BOFASA operates in two distinct phases:

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

A separate, lightweight subcommand, `bofasa domscape`, builds a Pfam domain by genome matrix without running `prep`/`run` (see [below](#standalone-pfam-domain-matrix-bofasa-domscape)).

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
- `-mm, --max-memory`: Optional cap, in GB, on virtual memory (address space, via RLIMIT_AS) - not resident RAM. It is inherited by every external tool bofasa runs and usually needs to be well above actual RAM use (default: no limit)

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
- Use PyHMMER to search for Pfam-A domains in proteins (using Pfam gathering thresholds)
- Select a set of mostly non-overlapping domains per protein:
   - Domain hits are processed from highest to lowest bitscore.
   - A hit is accepted unless it overlaps already accepted domains by ≥10% of its length or by ≥ `-ml` residues.
   - When an accepted hit has a minor (<10%) overlap with previously accepted domains, the overlapping residues are trimmed from its ends so that each accepted domain yields exactly one, correctly labelled, chunk. If trimming would leave nothing, or another accepted domain is nested within the hit, the hit is skipped.
- Create "chopped-up" protein FASTA per input genome where each protein is split up into "chunks" based on domain boundary coordinates: domain chunks plus the inter-domain regions between them (chunks shorter than `-ml` are not kept).

**Key Parameters:**
- `-ml, --min-length` (default: 20): Minimum length of a domain/inter-domain chunk
- `-sds, --skip-domain-splitting`: Skip domain-based splitting of proteins
- `-ntdo, --no-trim-domain-overlaps`: Do not trim minor overlaps between accepted domains (restores v1.2.2 behavior)

**Reasons:**
- Increases sensitivity by analyzing homology at domain resolution
- Accounts for domain shuffling/fusion/loss events
- Downstream we will be using OrthoFinder to determine coarse domain ortholog groups, which standardizes for differences in protein-chunk length

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
- `-gs, --genomad-splits` (default: 8): Number of chunks geNomad splits its database into (geNomad `--splits`); more splits lower memory usage at the cost of speed

**Effects:**
- `-rg` alone: Annotates MGEs but keeps them in original genomes
- `-rg` + `-emg`: Creates separate genome files for each phage/plasmid

**Rationale:**
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

**General Parameters (apply to all of `bofasa run`):**
- `-c, --threads` (default: 4): Number of threads/processes
- `-mm, --max-memory`: Optional cap, in GB, on virtual memory (not resident RAM); same behavior as in `bofasa prep` (default: no limit)

> [!NOTE]
> When MGEs are present (from `bofasa prep -emg`), they are **excluded** from OrthoFinder to prevent biasing the core ortholog group structure. MGEs are integrated in Step 2.

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
   - Filter hits: require E-value < 1e-3
   - Optionally (`-mmc`), also require a minimum query **and** subject coverage
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

**Key Parameters:**
- `-mmc, --mge-min-coverage` (default: off): Minimum query and subject coverage (%, 0-100) of MGE-vs-MGE DIAMOND hits used to cluster MGE protein chunks. Without a coverage filter, chained clusters spanning several OGs are left as singletons in step 7, which keeps false-positive assignments low.

## Step 3: Phylogenetic Refinement of Domain-Resolution Ortholog Groups

```
Coarse Domain-Resolution Ortholog Groups from OrthoFinder
  ↓
Per-Group Multiple Sequence Alignment using FAMSA (or MUSCLE super5)
  ↓
Phylogenetic Tree Construction using FastTree2
  ↓
Tree Analysis to Split Paralogs
  ↓
Refined Domain Ortholog Groups
```

**Algorithm:**
For each **coarse domain-resolution ortholog group**:
   - Build multiple sequence alignment (default: FAMSA; MUSCLE in `super5` mode can be requested via `--aligner muscle`).
   - Trim MSA using trimAl (in `strict` mode).
   - Construct phylogenetic tree (default: FastTree2).
   - Apply midpoint-rooting. If `-rs N` flag is specified, then root the tree at _N_-1 randomly chosen nodes and once using midpoint. The random nodes are drawn with a fixed seed per coarse ortholog group (derived from its name), so the rootings tried do not depend on the number of threads or the order in which jobs run. Each (tree, rooting) combination is processed as an independent job across `-c` processes, largest trees first.
   - Partition phylogeny into **refined domain-resolution ortholog groups** using bofasa's recursive algorithm.
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
   - Select the minimum-score subtrees greedily in preorder, skipping any subtree nested inside one already selected, so the selected subtrees do not overlap

#### Handle remaining sequences:
- Create "complement" partition (sequences not in best subtrees)
- ***If complement has ≥2 samples, recursively split it***
- Combine best subtrees + complement splits
  
### Details on additional refinement steps:

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

- Uses largest-to-smallest node traversal
- Calculates fixation index (FST) between split groups
- If FST < threshold → merge groups back together (default threshold = 0.1)
- Prevents over-splitting of true orthologs

**Key Parameters:**
- `-spr, --skip-phylo-refine`: Skip phylogenetic refinement entirely
  - Uses OrthoFinder results directly
  - Much faster but less accurate
  
- `-fic`, `--fixation-index-cutoff` (default: 0.1): Minimum FST to accept a split
  - Lower values (0.0-0.09): More re-merging
  - Higher values (0.11-1.0): Less re-merging
  
- `-smb, --skip-merge-back`: Skip merge-back assessment
  - Faster but may over-split ortholog groups
  
- `-rs, --rooting-seeds` (default: 1): Number of outgroup rooting attempts
  - Higher values: More robust rooting but slower
  - 1 (default): Fastest, uses midpoint rooting
  - 100: Tries midpoint rooting + 99 random nodes as roots to see if they result in higher quality refined ortholog groups
  - Best partitioning is selected based on scoring:
     - Sort by: (1) tree score, (2) number of groups, (3) total branch length and choose partition with lowest combined score
     - Ties are broken by the order rootings are evaluated: the random rootings first, midpoint rooting last
  
- `-al, --aligner`: Multiple sequence aligner to use: `famsa` (default) or `muscle` (MUSCLE super5). Note that MUSCLE super5 alignments are not deterministic (see https://github.com/rcedgar/muscle/issues/95), so runs using `muscle` are not exactly reproducible; FAMSA is deterministic.

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
1. For each protein, count how many chunks it has in each refined domain ortholog group (DOG)
2. For each pair of proteins sharing at least one DOG, calculate the Jaccard similarity of their DOG compositions:
   ```
   J(A,B) = |DOGs(A) ∩ DOGs(B)| / |DOGs(A) ∪ DOGs(B)|
   ```
   where _DOGs(A)_ is the multiset of DOGs of protein _A_, with each DOG counted once per chunk of _A_ in it. For multisets, ∩ keeps the smaller and ∪ the larger count of each DOG. For example, if _A_ has two chunks in DOG _X_ and _B_ has one, _X_ contributes 1 to the intersection and 2 to the union. DOGs that are largely inter-domain regions are excluded (see note below).
3. Keep pairs with J ≥ threshold. The threshold is `-dj` (default: 0.5), but is lowered to `-dj` × (1 − _c_) when the two proteins share a single-copy DOG (never more than one copy per genome) found in a proportion _c_ of genomes, using the most conserved such DOG. Proteins sharing a highly conserved single-copy DOG can therefore be grouped even if the rest of their domain composition differs.
4. Cluster the resulting protein similarity graph (edge weights = J) with MCL (inflation 1.2; independent of `-mi`) = coarse protein ortholog groups

### Part B: Phylogenetic Refinement of Multi-Copy Groups
For protein ortholog groups meeting refinement criteria (≥2 copies in any genome, ≥2 genomes, ≥4 proteins total):

1. **Domain-based distance calculation**:
   - Represent each protein as a vector of domain ortholog group counts
   - Calculate pairwise cosine distances between protein domain vectors
   - Create cosine-distance matrix in PHYLIP format

2. **Coarse protein ortholog group tree construction**:
   - Build phylogenetic tree using FastME

3. **Tree-based splitting**:
   - Root tree at midpoint
   - Apply recursive splitting at duplication nodes (uses same recursive function as for domain-resolution OG splitting)
       - Note, however, that there is a slight difference for domain-resolution and protein-resolution recursive splitting. Namely, for proteins, we stop only if the branch length sum = 0 (identical) whereas for domains we stop if branch length ≤ threshold. This is in order to account for pseudovalues added by FastTree 2.

4. **Monophyletic enforcement**:
   - Check if each split group forms a monophyletic clade
   - If non-monophyletic, further split into largest monophyletic sub-clades
   - Ensures evolutionary coherence within final ortholog groups

5. **Outlier artifact removal**:
   - For each protein in a split group:
     - Calculate max Jaccard similarity (of DOG multisets, as in Part A) to proteins within the group
     - Calculate max Jaccard similarity to proteins outside the group (elsewhere in the coarse group's tree)
     - If external similarity ≥ internal similarity → remove as singleton
   - Prevents misplaced proteins from remaining in inappropriate groups - could result from nesting pattern in FastME tree that arises akin to similar artifacts in neighbor-joining trees for instance.

> [!NOTE]
> When calculating Jaccard similarity between proteins based on domain-resolution ortholog groups (DOGs) for outlier detection, DOGs that are largely (>80%) inter-domain regions (IDR) are included to provide higher resolution. However, for initial aggregation of proteins into coarse protein ortholog groups (Part A), these largely-IDR DOGs are not counted as they might contribute noise.

**Key Parameters:**
- `-dj, --dog-jaccard` (default: 0.5): Jaccard similarity threshold for initial clustering (relaxed for proteins sharing conserved single-copy DOGs, see Part A)
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
- **Total Distinct Neighbors**: Total number of unique OGs that ever appear as neighbors to the focal OG across ***all*** contexts.

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
1. Align all sequences in each ortholog group (using the aligner set by `-al`; FAMSA by default)
2. Determine consensus sequence + profile HMM for each ortholog group

## Optional Step 8: Core Genome Alignment

### Premise:

Construct a concatenated multi-FASTA alignment of the strict or loose single-copy-core genome for downstream phylogenomics. It uses domain-resolution ortholog groups. 

### Key Parameter:
- `-cg, --core-genome`: Construct concatenated core genome alignment
- `-ns, --near-scc-prop` (default: 0.95): Minimum proportion of genomes for a domain ortholog group to be considered part of the single-copy-core.
  - Lower values (0.80-0.90): More permissive, more genes
  - Higher values (0.98-1.00): Stricter core, fewer genes

> [!IMPORTANT]
> An ortholog group can still be considered part of the loose single-copy-core (`-ns` <1.0) when genomes have multiple copies of it. These are ignored similar to the absence of the gene in other genomes. For instance, if your dataset has 100 genomes and 97 of the genomes have the focal ortholog group in single-copy, but 2 genomes lack it and 1 genome has two copies of the ortholog group, if `-ns` >= 0.97, then it is still considered part of the scc and treated as absent in that 1 genome with paralogs. **Also, if you are _not_ using a strict core genome (`-ns` set to 1.0), we recommend that you use a partition-based approach for phylogeny modeling and using the individual ortholog group alignments folder as input to IQ-TREE instead of the concatenated multi-FASTA file.**
  
### Algorithm:
1. Identify single-copy core ortholog groups (SCC-DOGs)
   - Present in ≥95% of genomes by default unless `-ns` modified.
2. Align each SCC-DOG (using the aligner set by `-al`; FAMSA by default)
3. Concatenate alignments and output multi-FASTA suitable for phylogenomics

---

# Database Setup (`bofasa setup`)

`bofasa setup` downloads the databases used by bofasa (needs to be run only once):

- **Pfam-A**: By default, the current Pfam release is downloaded. A specific release can be requested with `-pv, --pfam-version` (currently `37.4`, the release used to generate the test case's expected results); the download is validated against the expected number of HMMs for that release.
- **ISfinder** protein sequences: pinned to a fixed commit (October 2020) of the [ISfinder-sequences](https://github.com/thanhleviet/ISfinder-sequences) repository.
- **geNomad** database: only needed for `bofasa prep -rg`; can be skipped with `--skip-genomad`.

Pfam-A and ISfinder are downloaded with curl. Files of 64 MB or more are fetched over several parallel connections in byte-range segments (`-ds, --download-segments`, default: 8; `1` disables this) when the server supports range requests, and the reassembled file is checked against the size the server reports. The geNomad database is downloaded by `genomad download-database` and is not affected by `-ds`.

> [!NOTE]
> Domain splitting, and therefore all downstream results, depends on the Pfam release. Use the same Pfam release (e.g. via `-pv`) when comparing or reproducing runs.

---

# Standalone Pfam Domain Matrix (`bofasa domscape`)

A lightweight alternative to the full pipeline for quickly surveying domain content across many genomes. It is independent of `prep`/`run` and has no upper limit on the number of genomes.

```
Genomes (FASTA/GenBank), proteomes, or annotation dirs (Prokka/Bakta)
  ↓
Gene Calling (if FASTA genomes)
  ↓
PyHMMER search against Pfam-A (gathering thresholds), proteomes searched in batches
  ↓
Pfam domain by genome matrix (number of proteins carrying each domain)
```

**Algorithm:**
1. Gene calling (Pyrodigal/Prodigal) for genomes given in FASTA format; proteins are taken directly from GenBank files, annotation directories, or input proteomes.
2. Pfam HMMs are loaded once and proteomes are searched in batches (`-bs`).
3. Each protein is counted once per Pfam domain, even if it carries multiple copies of that domain.
4. Per-genome results are written to disk as they complete, so an interrupted run resumes when re-run with the same output directory.

**Key Outputs:**
- `Pfam_Protein_Count_Matrix.tsv`: Pfam domains (rows) by genomes (columns), with the number of proteins carrying each domain
- `Pfam_Domain_Info.tsv`: Per-domain name, number of genomes with the domain, and total proteins with the domain
- `Pfam_Domain_Proteins.tsv`: Listing of each protein carrying each domain in each genome

**Key Parameters:**
- `-i, --input-genomes` / `-p, --input-proteomes` / `-a, --annotation-dirs`: Inputs (at least one required)
- `-ml, --min-length` (default: 20): Minimum protein length kept from genomes/annotation directories (not applied to input proteomes)
- `-bs, --batch-size` (default: 100): Number of proteomes searched against Pfam together; larger is faster but uses more memory

---

# Reproducibility

As of v1.3.0, two runs on the same inputs, with the same Pfam release and the same options, are designed to give identical results:

- FAMSA (the default aligner) is deterministic; MUSCLE super5 (`-al muscle`) is not (see https://github.com/rcedgar/muscle/issues/95)
- Sequences are sorted before alignment, and protein pairs are sorted before MCL, so tool inputs do not depend on processing order
- Random rootings (`-rs` > 1) use a fixed seed per coarse ortholog group, so they do not depend on the number of threads or job order

---

## Key Citations

When using BOFASA, please cite:

> Salamzade, R., Yaffe, E., Kottapalli, A., Kalan, L., & Relman, D. (2026). bofasa enables high-quality orthology inference across diverse multi-species bacterial genome sets. bioRxiv.

And the underlying tools:
- **OrthoFinder**: Emms, D.M. & Kelly, S. (2019) Genome Biology
- **geNomad**: Camargo, A.P. et al. (2023) Nature Biotechnology
- **FastTree**: Price, M.N. et al. (2010) PLoS One (domain-level phylogeny)
- **FastME**: Lefort, V. et al. (2015) Molecular Biology and Evolution (protein-level phylogeny)
- **FAMSA**: Deorowicz, S. et al. (2016) Scientific Reports
- **MUSCLE**: Edgar, R.C. (2022) Nature Communications
