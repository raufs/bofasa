# bofasa

**B**acterial **O**rthology **F**inding **A**nd **S**yntenic **A**nalysis (**bofasa**)

<p align="center">
<img src="https://github.com/user-attachments/assets/11fa592f-656e-4578-93b8-17f8dc36f6c2" width="400">
</p>

bofasa is specifically designed for investigating orthology between multiple-species of bacteria. It prioritizes high-quality orthology inference at the expense of throughput, designed to run on between 4 and 200 genomes.

For single species analyses, we recommend [Panaroo](https://github.com/gtonkinhill/panaroo) or [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN), which for such cases, offers considerable advantages in terms of speed, accuracy, and scalability.

Comparable and good alternatives to bofasa include [PIRATE](https://github.com/SionBayliss/PIRATE) and [SCARAP](https://github.com/SWittouck/SCARAP). PIRATE performs similarly to bofasa when applied to genomes representative of different species from a single genus when using default parameters. SCARAP behaves more similar to non-bacteria specific multi-species orthology inference software such as OrthoFinder and SonicParanoid and is really fast and could be well-suited for more comprehensive and large-scale analysis. 

## Quick Start

bofasa provides a unified command-line interface with two main subcommands:

### 1. Prepare Input Data
```bash
bofasa prep -i genome1.fasta genome2.fasta genome3.fasta genome4.fasta -o prep_output
```

### 2. Run Analysis
```bash
bofasa run -i prep_output -o analysis_output -c 8
```

## Installation

### Bioconda

**_Coming soon ..._**

### Conda 

```bash
# git clone repository
git clone https://github.com/your-repo/bofasa.git
cd bofasa

# create conda environment and activate it
conda env create -f bofasa_env.yml -n bofasa
conda activate bofasa

# pip install the program
pip install -e .
```

## Usage

### Main Interface

bofasa provides a single entry point with intuitive subcommands:

```
██████╗  ██████╗ ███████╗ █████╗ ███████╗ █████╗
██╔══██╗██╔═══██╗██╔════╝██╔══██╗██╔════╝██╔══██╗
██████╔╝██║   ██║█████╗  ███████║███████╗███████║
██╔══██╗██║   ██║██╔══╝  ██╔══██║╚════██║██╔══██║
██████╔╝╚██████╔╝██║     ██║  ██║███████║██║  ██║
╚═════╝  ╚═════╝ ╚═╝     ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝

BOFASA: Bacterial Ortholog Finder and Synteny Analysis

A comprehensive tool for identifying ortholog groups and analyzing
syntenic conservation in bacterial genomes representing multiple
species.

Commands:
  setup    Set up annotation databases
  prep     Prepare genomic data for analysis
  run      Run the main BOFASA analysis pipeline

For detailed help on any command, use: bofasa <command> --help

Examples:
  bofasa prep -i genome1.fna genome2.fna -o prep_output/
  bofasa run -i prep_output/ -o analysis_results/ -c 8

For more information, visit: https://github.com/raufs/bofasa
```

> [!NOTE]
> The concepts in bofasa were developed over many years and the initial release (v1.1.0) developed without the use of AI. For version v1.2.0, we used AI to largely restructure the code and in the process caught some bugs and added a couple optimizations/bells and whistles. Moving forward, we will likely use some AI because its engraved into search engines but will largely try to refrain from relying on it to have a sense of the code so that it is not yet another black box.

### Subcommands

#### `bofasa prep` - Prepare Input Data
Processes input genomes (FASTA or GenBank files) for analysis.

```bash
bofasa prep -i <genome1.fasta> <genome2.gbk> ... -o <output-dir> [OPTIONS]
```

**Required Arguments:**
- `-i, --input-genomes`: Input genomes (FASTA or GenBank files)
- `-o, --output-dir`: Output directory for prepared data

**Optional Arguments:**
- `-c, --threads`: Number of threads (default: 4)
- `-gcm, --gene-calling-method`: Gene calling method (pyrodigal/prodigal, default: pyrodigal)
- `-l, --locus-tag-length`: Length of locus tags (default: 3)
- `-m, --meta-mode`: Use meta-mode for gene calling
- `-rlt, --rename-locus-tags`: Rename locus tags in GenBank files
- `-rg, --run-genomad`: Run genomad for phage/plasmid annotation
- `-emg, --extract-mge-genomes`: Extract mobile genetic element (phage/plasmid) genomes identified by genomad (requires -rg)
- `-mm, --max-memory`: Memory limit in GB (default: 32)

#### `bofasa run` - Execute Analysis
Performs the main orthology analysis on prepared data.

```bash
bofasa run -i <bofasa-prep-dir> -o <output-dir> [OPTIONS]
```

**Required Arguments:**
- `-i, --bofasa-prep-dir`: Input directory produced by bofasa prep
- `-o, --output-dir`: Output directory

**Optional Arguments:**
- `-sr, --surrounding-bp`: Base pairs for syntenic analysis (default: 10000)
- `-ogc, --og-consensus`: Determine consensus sequences
- `-cg, --core-genome`: Construct core genome alignment
- `-dj, --dog-jaccard`: Jaccard index threshold (default: 0.25)
- `-fic, --fixation-index-cutoff`: Fixation index cutoff (default: 0.25)
- `-smb, --skip-merge-back`: Skip merge back assessment
- `-spr, --skip-phylo-refine`: Skip phylogenetic refinement
- `-rs, --rooting-seeds`: Rooting seeds (default: 1)
- `-qa, --quality-alignments`: Use high-quality alignments
- `-us, --ultra-sens`: Use ultra-sensitivity mode
- `-mi, --mcl-inflation`: MCL inflation parameter (default: 1.2)
- `-ns, --near-scc-prop`: Near SCC proportion (default: 0.95)
- `-c, --threads`: Number of threads (default: 4)
- `-mrd, --max-recursion-depth`: Max recursion depth (default: 5000)
- `-mm, --max-memory`: Memory limit in GB (default: 32)

### Examples

```bash
# Prepare input data
bofasa prep -i genome1.fasta genome2.gbk genome3.fasta -o prep_output -c 8

# Run analysis with custom parameters
bofasa run -i prep_output -o analysis_output -c 8 -mm 64 -cg -ogc

# Get help
bofasa --help
bofasa prep --help
bofasa run --help

# Get version
bofasa --version
```

## Additional Scripts

bofasa also provides several additional utility scripts:

- `extract_og_pairs` - Extract ortholog group pairs for method comparison
- `extract_scc_og_pairs` - Extract single-copy-core ortholog group pairs
- `compare_og_pairs` - Compare ortholog group pairs between methods
- `determine_og_freqs` - Determine ortholog group frequencies
- `print_og_itol_matrix` - Print ortholog group matrices for iTOL visualization

**Note:** GenBank processing and Prodigal gene calling functionality is now integrated into the main `bofasa prep` command and no longer requires separate scripts.

## Documentation

For detailed information on bofasa, see the `ALGORITHM.md` document.

## Citation

> ***High-quality and automated inference of ortholog groups across multiple bacterial species using bofasa.***
> Rauf Salamzade, Eitan Yaffe, Aamuktha Kottapalli, David Relman, Lindsay Kalan, 2026.

Please also consider citing both ***OrthoFinder 2*** and ***geNomad*** which are used for determination of coarse domain ortholog groups and the annotation of phages/plasmids, respectively.

## License

```
BSD 3-Clause License

Copyright (c) 2024, Rauf Salamzade

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```
