# BOFASA Algorithm Documentation

## Table of Contents
1. [Overview](#overview)
2. [Workflow Architecture](#workflow-architecture)
3. [Preparation Phase (`bofasa prep`)](#preparation-phase)
4. [Analysis Phase (`bofasa run`)](#analysis-phase)
5. [Algorithm Details](#algorithm-details)
6. [Parameter Effects](#parameter-effects)
7. [Best Practices](#best-practices)

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

## Preparation Phase (`bofasa prep`)

The preparation phase processes raw genome data into a format suitable for orthology analysis.

### Step-by-Step Process

#### 1. Genome Input Processing

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

#### 2. Protein Extraction
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

#### 3. Domain Identification (via PyHMMER+Pfam)
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

#### 4. Optional: Mobile Genetic Element Detection

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

## Analysis Phase (`bofasa run`)

The analysis phase performs hierarchical ortholog inference using the prepared domain and protein data.

### Step 1: Inference of Coarse Domain-Resolution Ortholog Groups using OrthoFinder

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

### Step 2: MGE Integration (_auxiliary_; if `-emg` used in prep)

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

### Step 2: Phylogenetic Refinement of Domain-Resolution Ortholog Groups

```
Coarse Domain-Resolution Ortholog Groups from OrthoFinder
  ↓
Per-Group Multiple Sequence Alignment
  ↓
Phylogenetic Tree Construction
  ↓
Tree Analysis & Splitting
  ↓
Refined Domain Ortholog Groups
```

**Algorithm:**
1. For each coarse ortholog group:
   - Build multiple sequence alignment (MUSCLE)
   - Construct phylogenetic tree (FastTree2 or IQ-TREE)
   - Root the tree using multiple outgroup strategies
   - Identify duplication nodes (gene duplications)
   - Split groups at duplication nodes
   - Assign sequences to refined ortholog groups

2. Fixation index calculation:
   - For each potential split, calculate fixation index (Fst)
   - Fst measures genetic differentiation between groups
   - Higher Fst = stronger evidence for distinct ortholog groups

3. Merge-back assessment (optional):
   - Check if split groups should be reunited
   - Prevents over-splitting of true orthologs

**Key Parameters:**
- `-spr, --skip-phylo-refine`: Skip phylogenetic refinement entirely
  - Uses OrthoFinder results directly
  - Much faster but less accurate
  
- `-fic, --fixation-index-cutoff` (default: 0.25): Minimum Fst to accept a split
  - Lower values (0.1-0.2): More aggressive splitting
  - Higher values (0.3-0.5): Conservative splitting
  
- `-smb, --skip-merge-back`: Skip merge-back assessment
  - Faster but may over-split ortholog groups
  
- `-rs, --rooting-seeds` (default: 1): Number of outgroup rooting attempts
  - Higher values: More robust rooting but slower
  - 1: Fastest, usually sufficient
  - 3-5: Better for complex paralog situations
  
- `-qa, --quality-alignments`: Use high-quality alignment settings
  - MAFFT instead of MUSCLE
  - IQ-TREE instead of FastTree2
  - Much slower but better for challenging cases

**Effects:**
- **Skip phylogenetic refinement** (`-spr`):
  - ✓ 5-10x faster
  - ✓ Deterministic results
  - ✗ May include paralogs in ortholog groups
  - ✗ Lower accuracy for recently duplicated genes
  - **Use when**: Species are very closely related OR speed is critical

- **Lower fixation index cutoff** (e.g., 0.15):
  - ✓ More aggressive paralog detection
  - ✓ Better separation of gene families
  - ✗ May split true orthologs
  - ✗ More fragmentation
  - **Use when**: Analyzing genomes with many recent duplications

- **Higher fixation index cutoff** (e.g., 0.35):
  - ✓ Conservative splitting
  - ✓ Fewer false splits
  - ✗ May retain some paralogs
  - **Use when**: Analyzing ancient orthologs with strong conservation

- **Skip merge-back** (`-smb`):
  - ✓ Faster
  - ✗ May over-split reciprocal best hit relationships
  - **Use when**: Groups are well-separated OR computational time is limited

- **Multiple rooting seeds** (e.g., `-rs 5`):
  - ✓ More robust to outgroup selection
  - ✓ Better handling of complex gene families
  - ✗ Significantly slower (linear with number of seeds)
  - **Use when**: Working with complex multi-domain proteins OR paralog-rich families

- **Quality alignments** (`-qa`):
  - ✓ Better phylogenetic accuracy
  - ✓ More reliable for divergent sequences
  - ✗ 10-50x slower depending on group size
  - **Use when**: Final publication-quality analysis OR very divergent species

### Step 3: Protein Ortholog Group Determination

```
Refined Domain Ortholog Groups
  ↓
Protein-to-Domain Mapping
  ↓
Domain Jaccard Similarity Calculation
  ↓
Graph-based Protein Clustering
  ↓
Coarse Protein Ortholog Groups
  ↓
Phylogenetic Refinement (for multi-copy groups)
  ↓
Final Protein Ortholog Groups (POGs)
```

**Algorithm:**

#### Part A: Initial Jaccard-Based Clustering
1. For each protein, identify its domain ortholog group memberships
2. Calculate Jaccard similarity between proteins:
   ```
   J(A,B) = |DOGs(A) ∩ DOGs(B)| / |DOGs(A) ∪ DOGs(B)|
   ```
3. Create protein similarity graph (edges = similarity ≥ threshold)
4. Identify connected components = coarse protein ortholog groups
5. Handle single-domain proteins and proteins without domains

#### Part B: Phylogenetic Refinement of Multi-Copy Groups
For protein ortholog groups meeting refinement criteria (≥2 copies in any genome, ≥2 genomes, ≥4 proteins total):

1. **Domain-based distance calculation**:
   - Represent each protein as a vector of domain ortholog group counts
   - Calculate pairwise cosine distances between protein domain vectors
   - Create distance matrix in PHYLIP format

2. **Neighbor-joining tree construction**:
   - Build phylogenetic tree using FastME
   - Tree reflects evolutionary relationships based on domain architecture

3. **Tree-based splitting**:
   - Root tree at midpoint
   - Apply recursive splitting at duplication nodes (same as domain OG refinement)
   - Identify paralogs vs orthologs based on tree topology

4. **Outlier artifact removal**:
   - For each protein in a split group:
     - Calculate max Jaccard similarity to proteins within the group
     - Calculate max Jaccard similarity to proteins outside the group
     - If external similarity ≥ internal similarity → remove as singleton
   - Prevents misplaced proteins from remaining in inappropriate groups

5. **Monophyletic enforcement**:
   - Check if each split group forms a monophyletic clade
   - If non-monophyletic, further split into largest monophyletic sub-clades
   - Ensures evolutionary coherence within final ortholog groups

**Key Parameters:**
- `-dj, --dog-jaccard` (default: 0.25): Jaccard similarity threshold for initial clustering
  - Lower values (0.1-0.2): More permissive, allows domain rearrangements
  - Higher values (0.3-0.5): Stricter, requires more domain conservation

**Effects:**
- **Lower Jaccard threshold** (e.g., 0.15):
  - ✓ Groups proteins with domain rearrangements
  - ✓ Better for multi-domain protein families
  - ✓ Captures domain shuffling events
  - ✗ May group functionally distinct proteins
  - ✗ Requires more refinement splitting
  - **Use when**: Studying multi-domain proteins OR domain evolution

- **Higher Jaccard threshold** (e.g., 0.4):
  - ✓ Stricter domain architecture conservation
  - ✓ More functionally coherent groups
  - ✓ Less refinement needed
  - ✗ May split orthologs with minor domain differences
  - ✗ More singleton proteins
  - **Use when**: Analyzing single-domain proteins OR requiring high confidence

**Why Phylogenetic Refinement is Necessary:**

Even after Jaccard-based clustering, protein ortholog groups can contain paralogs because:
1. Paralogs may share similar domain architectures (especially recent duplications)
2. Jaccard threshold may be permissive to avoid splitting true orthologs
3. Domain rearrangements can create misleading similarity patterns

The phylogenetic refinement step uses evolutionary relationships to distinguish:
- **True orthologs**: Descended from speciation events, form species-congruent clades
- **Paralogs**: Descended from duplication events, show within-species clustering

### Step 4: Syntenic Context Analysis

```
Protein Ortholog Groups
  ↓
Extract Genomic Neighborhoods
  ↓
Calculate Neighborhood Conservation
  ↓
Compute Context Entropy
  ↓
Ortholog Groups + Synteny Metrics
```

**Algorithm:**
1. For each ortholog group occurrence:
   - Extract upstream and downstream neighbors
   - Identify their ortholog group assignments
   - Calculate neighborhood conservation metrics

2. Metrics calculated:
   - **Neighborhood entropy**: Shannon entropy of neighboring OG composition
     - Low entropy = conserved synteny
     - High entropy = variable synteny
   - **Neighbor presence/absence matrix**: Which OGs appear as neighbors
   - **Directional conservation**: Strand orientation patterns

**Key Parameters:**
- `-sr, --surrounding-bp` (default: 10000): Base pairs to analyze around each gene
  - Smaller values (5000-7500): Immediate neighbors only
  - Larger values (15000-25000): Extended chromosomal context

**Effects:**
- **Smaller surrounding region** (e.g., 5000 bp):
  - ✓ Faster computation
  - ✓ Focus on immediate operonic structure
  - ✗ May miss larger syntenic patterns
  - ✗ Less context for sparse genomes
  - **Use when**: Analyzing operon structure OR gene-dense genomes

- **Larger surrounding region** (e.g., 20000 bp):
  - ✓ Captures broader chromosomal organization
  - ✓ Better for sparse genomes
  - ✓ More robust to local rearrangements
  - ✗ Slower computation
  - ✗ May include unrelated genes
  - **Use when**: Studying genomic islands OR large-scale synteny

### Step 5: Final Report Generation

```
All Analysis Results
  ↓
Data Integration
  ↓
Excel Spreadsheet Generation
```

**Outputs:**
- Multi-sheet Excel workbook with:
  - Ortholog group summary (presence/absence, copy numbers)
  - Synteny conservation metrics
  - Per-sample statistics
  - Functional annotations (if available)

### Step 6: Visualization

```
Synteny + Conservation Data
  ↓
Interactive Plotting
  ↓
HTML Visualization
```

**Output:**
- Interactive HTML plot showing:
  - X-axis: Ortholog group conservation (% genomes)
  - Y-axis: Neighborhood entropy (syntenic conservation)
  - Hover: Detailed OG information

### Optional Step 7: Consensus Sequence Generation

**Key Parameter:**
- `-ogc, --og-consensus`: Generate consensus sequences for each ortholog group

**Algorithm:**
1. Align all sequences in each ortholog group
2. Call consensus at each position (majority rule)
3. Create consensus FASTA file
4. Build profile HMMs for each ortholog group

**Effects:**
- ✓ Enables downstream searching/annotation
- ✓ Profile HMMs useful for finding orthologs in new genomes
- ✗ Significantly increases runtime (proportional to # of OGs)
- ✗ Large disk space requirement
- **Use when**: Planning to annotate new genomes OR need representative sequences

### Optional Step 8: Core Genome Alignment

**Key Parameter:**
- `-cg, --core-genome`: Construct concatenated core genome alignment

**Algorithm:**
1. Identify single-copy core ortholog groups (SCC-OGs)
   - Present in ≥95% of genomes (default)
   - Exactly one copy per genome
2. Align each SCC-OG
3. Concatenate alignments
4. Output multi-FASTA suitable for phylogenomics

**Key Parameters:**
- `-ns, --near-scc-prop` (default: 0.95): Minimum proportion of genomes for "core"
  - Lower values (0.80-0.90): More permissive, more genes
  - Higher values (0.98-1.00): Strict core, fewer genes

**Effects:**
- **Lower near-SCC proportion** (e.g., 0.85):
  - ✓ More genes in core genome alignment
  - ✓ More phylogenetic signal
  - ✗ May include non-universal genes

- **Higher near-SCC proportion** (e.g., 0.98):
  - ✓ Very strict core genome
  - ✓ Suitable for diverse species sets
  - ✗ Fewer genes, less phylogenetic signal

---

## Algorithm Details

### Domain Orthology vs Protein Orthology

BOFASA's hierarchical approach offers advantages over protein-only methods:

```
Traditional Approach:
Protein A1 [Domain1-Domain2-Domain3]  ─┐
Protein B1 [Domain1-Domain2]           ├─ May or may not cluster together
Protein C1 [Domain1-Domain3]          ─┘

BOFASA Approach:
Step 1 - Domain level orthology:
  Domain1: A1, B1, C1 → Ortholog Group 1
  Domain2: A1, B1     → Ortholog Group 2
  Domain3: A1, C1     → Ortholog Group 3

Step 2 - Jaccard-based protein clustering:
  Jaccard(A1, B1) = 2/3 = 0.67 → Link
  Jaccard(A1, C1) = 2/3 = 0.67 → Link
  Jaccard(B1, C1) = 1/3 = 0.33 → Link (if threshold ≤0.33)

Step 3 - Phylogenetic refinement (if multi-copy):
  Build tree based on domain composition distances
  Split at duplication nodes
  Remove misplaced proteins
  Enforce monophyletic property
```

**Advantages:**
- Handles domain shuffling elegantly
- More sensitive for multi-domain proteins
- Reduces false negatives from domain rearrangements
- Phylogenetic refinement catches paralogs missed by Jaccard clustering
- Two-stage approach balances sensitivity and specificity

**Trade-offs:**
- More complex pipeline
- Longer runtime
- Depends on domain prediction accuracy

### Phylogenetic Refinement Details

The phylogenetic refinement step is critical for distinguishing orthologs from paralogs at both the domain and protein levels.

#### Domain-Level Phylogenetic Refinement

Applied to coarse domain ortholog groups from OrthoFinder:

**Rooting Strategy:**
```
Unrooted Tree
  ↓
Try multiple rooting strategies:
  1. Minimize tree depth variance
  2. Balance by species representation
  3. Outgroup selection (if available)
  ↓
Select best rooting (highest confidence)
```

**Duplication Detection:**
```
For each internal node:
  1. Check if children have overlapping species
  2. If yes → duplication node
  3. Calculate Fst between child clades
  4. If Fst ≥ threshold → accept split
```

**Fixation Index (Fst) Interpretation:**
- **Fst = 0**: No differentiation (likely one ortholog group)
- **Fst = 0.1-0.25**: Moderate differentiation (borderline)
- **Fst = 0.25-0.5**: Strong differentiation (likely distinct orthologs)
- **Fst > 0.5**: Very strong differentiation (clear paralogs)

#### Protein-Level Phylogenetic Refinement

Applied to coarse protein ortholog groups from Jaccard clustering:

**Refinement Criteria:**
Protein ortholog groups are refined if they meet ALL of:
- At least 2 copies present in any single genome (indicates potential paralogs)
- Present in at least 2 different genomes
- Contains at least 4 proteins total

**Distance Calculation:**
Unlike domain-level refinement (which uses sequence alignment), protein-level refinement uses domain architecture:
```
For proteins P1 and P2:
  1. Create domain composition vectors:
     V(P1) = [count(DOG_1), count(DOG_2), ..., count(DOG_n)]
     V(P2) = [count(DOG_1), count(DOG_2), ..., count(DOG_n)]
  
  2. Calculate cosine distance:
     dist(P1, P2) = 1 - (V(P1) · V(P2)) / (||V(P1)|| × ||V(P2)||)
```

This approach is faster than full sequence alignment and captures evolutionary relationships based on domain gain/loss events.

**Tree Construction and Splitting:**
```
Distance Matrix (PHYLIP format)
  ↓
FastME (neighbor-joining)
  ↓
Midpoint rooting
  ↓
Recursive splitting at duplication nodes
  ↓
Outlier artifact removal
  ↓
Monophyly enforcement
  ↓
Final refined protein ortholog groups
```

**Outlier Artifact Removal:**

Prevents tree artifacts from creating spurious groups:
```python
For each protein P in split group G:
  max_internal_similarity = max(Jaccard(P, X) for X in G, X ≠ P)
  max_external_similarity = max(Jaccard(P, Y) for Y not in G)
  
  if max_external_similarity >= max_internal_similarity:
    # P is more similar to proteins outside G than inside
    # Likely misplaced due to tree artifact → extract as singleton
    remove P from G
```

This catches cases where midpoint rooting or tree topology artifacts cause unrelated proteins to cluster together.

**Monophyly Enforcement:**

Ensures each final ortholog group forms a coherent evolutionary unit:
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

### Syntenic Context Entropy

Neighborhood entropy quantifies syntenic conservation:

```
H = -Σ(p_i * log(p_i))

where p_i = proportion of times OG_i appears as a neighbor
```

**Interpretation:**
- **H ≈ 0**: Highly conserved synteny (same neighbors always)
- **H = 1-2**: Moderate conservation (few common arrangements)
- **H > 3**: High variability (many different neighbors)

**Biological Interpretation:**
- Low entropy: Core metabolic genes, essential operons
- High entropy: Accessory genes, mobile elements, recently acquired

---

## Parameter Effects

### Computational Resource Parameters

#### Max Recursion Depth (`-mrd`)
- **Effect**: Limits phylogenetic tree recursion depth
- **When to adjust**: 
  - Increase (>5000) for very large ortholog groups
  - Decrease (<3000) if experiencing stack overflow

### Quality vs Speed Trade-offs

#### Fast Mode
```bash
bofasa run -i prep/ -o out/ -c 16 -spr
```
- Skips phylogenetic refinement
- ~5-10x faster
- Lower accuracy for paralogs
- **Use for**: Exploratory analyses, very close species

#### Balanced Mode (Default)
```bash
bofasa run -i prep/ -o out/ -c 16
```
- Standard phylogenetic refinement
- Good accuracy/speed balance
- **Use for**: Most analyses

#### High Quality Mode
```bash
bofasa run -i prep/ -o out/ -c 16 -qa -rs 5 -fic 0.30
```
- Best accuracy
- 10-50x slower
- **Use for**: Publication analyses, divergent species

---

## Best Practices

### Genome Selection

1. **Minimum genomes**: 4 (algorithm requirement)
2. **Maximum genomes**: 200 (performance degradation above this)
3. **Optimal range**: 10-50 genomes
4. **Species diversity**: 
   - Works best within a genus
   - Can handle multiple genera if relatively related
   - Not recommended for cross-phylum analyses

### Parameter Selection Guidelines

#### For Closely Related Species (ANI > 95%)
```bash
# Prep
bofasa prep -i *.fasta -o prep/

# Run
bofasa run -i prep/ -o out/ \
  -mi 1.4 \           # Higher MCL inflation
  -fic 0.20 \         # Lower fixation cutoff
  -dj 0.30 \          # Higher Jaccard threshold
  -sr 15000           # Larger syntenic window
```

**Note:** Closely related species often have recent gene duplications with high sequence similarity. The protein-level phylogenetic refinement step is particularly important here, as paralogs may share nearly identical domain architectures. The automatic refinement will split these based on phylogenetic relationships.

#### For Divergent Species (ANI < 85%)
```bash
# Prep  
bofasa prep -i *.fasta -o prep/

# Run
bofasa run -i prep/ -o out/ \
  -mi 1.1 \           # Lower MCL inflation
  -us \               # Ultra-sensitive DIAMOND
  -fic 0.30 \         # Higher fixation cutoff
  -dj 0.20 \          # Lower Jaccard threshold
  -qa \               # Quality alignments
  -rs 3               # Multiple rooting attempts
```

#### For Draft/Incomplete Genomes
```bash
# Prep
bofasa prep -i *.fasta -o prep/ -m  # Meta-mode

# Run
bofasa run -i prep/ -o out/ \
  -fic 0.25 \         # Standard fixation cutoff
  -sr 7500            # Smaller syntenic window
```

#### For Mobile Genetic Element Analysis
```bash
# Prep - Extract MGEs as separate entities
bofasa prep -i *.fasta -o prep/ -rg -emg

# Run - MGEs integrated via Step 1b
bofasa run -i prep/ -o out/ \
  -mi 1.0 \           # Lower inflation for distant MGE homologs
  -fic 0.15 \         # More aggressive splitting
  -sr 5000            # Smaller window (MGEs often lack synteny)
```

**MGE Analysis Notes:**
- MGEs are excluded from initial OrthoFinder run (Step 1)
- MGEs integrated in Step 1b via DIAMOND alignment against bacterial OGs
- Conservative approach: conflicting assignments → singletons
- Modified ortholog tables include MGE columns (`.ccds` suffix)
- Unassigned MGE proteins remain as singletons in modified tables

### Troubleshooting

#### Issue: OrthoFinder fails or hangs
**Causes:**
- Symlink issues on some filesystems
- Too many/too large input files
- Insufficient disk space

**Solutions:**
- BOFASA now copies files instead of symlinking (automatic)
- Temporary directories created: `OrthoFinder_Input_Bacterial_Genomes/`
- Check disk space in output directory
- Reduce number of input genomes if memory-limited

#### Issue: Too many singleton ortholog groups
**Causes:**
- MCL inflation too high
- Jaccard threshold too high
- Species too divergent

**Solutions:**
- Decrease `-mi` to 1.0-1.1
- Decrease `-dj` to 0.15-0.20
- Enable `-us` for ultra-sensitive mode
- Skip phylogenetic refinement (`-spr`) as a test

#### Issue: Ortholog groups contain obvious paralogs
**Causes:**
- MCL inflation too low (domain-level clustering too coarse)
- Fixation index cutoff too high (domain-level refinement too conservative)
- Phylogenetic refinement skipped (domain-level)
- Jaccard threshold too low (protein-level clustering too permissive)
- Multi-copy groups not meeting refinement criteria (<4 proteins or <2 genomes)

**Solutions:**
- Increase `-mi` to 1.3-1.5 (stricter domain clustering)
- Decrease `-fic` to 0.15-0.20 (more aggressive domain splitting)
- Add `-qa` for better tree quality (domain-level)
- Increase `-rs` to 3-5 (more robust domain-level rooting)
- Increase `-dj` to 0.30-0.35 (stricter protein clustering)

**Note:** Protein ortholog group refinement is automatic for qualifying groups (≥2 copies in any genome, ≥2 genomes, ≥4 total proteins). If paralogs persist, they may be in groups that don't meet these criteria or have very similar domain architectures.

#### Issue: Analysis too slow
**Causes:**
- Too many genomes
- Quality alignment mode enabled
- Large proteins/many domains

**Solutions:**
- Use `-spr` to skip phylogenetic refinement
- Remove `-qa` if enabled
- Reduce `-rs` to 1
- Consider subsampling genomes

#### Issue: High memory usage
**Causes:**
- Too many threads
- Large ortholog groups
- Many genomes

**Solutions:**
- Reduce `-c` threads
- Increase `-mm` memory limit (if available)
- Use `-spr` to reduce memory in refinement step

---

## Output Interpretation

### Key Output Files

1. **`Final_Results/Protein_Ortholog_Groups.tsv`**
   - Tab-delimited ortholog group matrix
   - Rows = ortholog groups
   - Columns = genomes
   - Values = comma-separated protein IDs

2. **`Final_Results/Orthogroup_Overview.xlsx`**
   - Multi-sheet Excel workbook
   - Sheet 1: Summary statistics per OG
   - Sheet 2: Presence/absence matrix
   - Sheet 3: Synteny metrics

3. **`Final_Results/Orthogroup_Conservation_vs_ContextEntropy.html`**
   - Interactive scatter plot
   - Identify core vs accessory genes
   - Find syntenic vs mobile genes

### Interpreting Results

#### Core vs Accessory Genes
- **Core genes** (present in >95% genomes):
  - Low conservation score variation
  - Usually low syntenic entropy
  - Essential functions

- **Accessory genes** (present in <50% genomes):
  - Variable presence
  - Often high syntenic entropy
  - Niche-specific adaptations

#### Syntenic Conservation Patterns
- **High conservation (low entropy)**:
  - Operonic organization
  - Co-regulated genes
  - Essential pathways

- **Low conservation (high entropy)**:
  - Recently acquired genes
  - Mobile genetic elements
  - Genes under relaxed selection

---

## Technical Implementation Notes

### File Organization and Workflow Optimizations

#### Unified Directory Structure
BOFASA uses a unified directory structure for all domain and coordinate files:
- **Before**: Separate subdirectories (`Bacterial_Genomes/`, `MGEs/`)
- **After**: Common directory with metadata-based filtering

**Benefits:**
- ✓ Simpler file organization
- ✓ Easier checkpoint recovery
- ✓ More maintainable code
- ✓ Logical separation maintained via `Sample_MGE_Metadata.txt`

#### File Handling for OrthoFinder
- **Approach**: Copies files instead of symlinks
- **Reason**: Better compatibility across filesystems
- **Implementation**: `shutil.copy2()` preserves metadata
- **Location**: Temporary directory `OrthoFinder_Input_Bacterial_Genomes/`

#### Large File Processing
**DIAMOND Output Sorting:**
- **Challenge**: DIAMOND output can be millions of lines
- **Solution**: Unix `sort` with multi-threading
- **Command**: `sort -t $'\t' -k12,12 -n -r --parallel={threads}`
- **Benefits**:
  - Disk-based external sorting (handles files larger than RAM)
  - Parallel processing for speed
  - Sorted by bitscore (descending) for quality-first assignment

**slclust Integration:**
- **Challenge**: Shell redirection operators (`<`, `>`) don't work with `subprocess.run()` without `shell=True`
- **Solution**: Explicit file handle redirection
- **Implementation**: `stdin=file_handle`, `stdout=file_handle`
- **Benefits**: Safer than shell=True, portable, proper error handling

#### MGE Integration Algorithm
**Homology Clustering:**
1. Reflexive DIAMOND alignment (MGE vs MGE)
2. Coverage filtering (≥50% query and subject coverage)
3. Single-linkage clustering (slclust)
4. Cluster-aware OG assignment

**Conservative Assignment Strategy:**
- Single OG consensus → Assign all cluster members
- Multiple OG assignments → Leave as singletons
- Rationale: Prioritizes accuracy over completeness

**Performance Optimizations:**
- Concatenated FASTA files for batch processing
- Separate DIAMOND databases for MGEs and bacterial genomes
- Bitscore-sorted results for quality-first processing
- Efficient set operations for sample filtering

### Code Quality and Error Handling

**Exception Handling Pattern:**
```python
try:
    # Operation
    pass
except Exception as e:
    log_object.error(f"Context: {str(e)}")
    log_object.error(traceback.format_exc())
    raise  # Re-raise with context logged
```

**Benefits:**
- Logs error context before propagating
- Preserves original exception and traceback
- Follows Python best practices (PEP8-compliant)

**OG Name Generation:**
- Centralized function: `generate_og_name(i)`
- Zero-padded formatting: `OG00001`, `OG00042`, etc.
- Ensures proper alphabetical sorting
- Used consistently across codebase

---

## Citation

When using BOFASA, please cite:

> Salamzade, R., Kottapalli, A., & Kalan, L. (2025). bofasa: high-quality orthology inference across multiple bacterial species.

And the underlying tools:
- **OrthoFinder**: Emms, D.M. & Kelly, S. (2019) Genome Biology
- **geNomad**: Camargo, A.P. et al. (2023) Nature Biotechnology
- **FastTree**: Price, M.N. et al. (2010) PLoS One (domain-level phylogeny)
- **FastME**: Lefort, V. et al. (2015) Molecular Biology and Evolution (protein-level phylogeny)
- **MUSCLE**: Edgar, R.C. (2022) Nature Communications

