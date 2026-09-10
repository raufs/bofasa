# bofasa
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bofasa/badges/version.svg)](https://anaconda.org/bioconda/bofasa)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bofasa/badges/latest_release_date.svg)](https://anaconda.org/bioconda/bofasa)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bofasa/badges/latest_release_relative_date.svg)](https://anaconda.org/bioconda/bofasa)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bofasa/badges/license.svg)](https://anaconda.org/bioconda/bofasa)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/bofasa/badges/downloads.svg)](https://anaconda.org/bioconda/bofasa)
[![Docker](https://img.shields.io/badge/Docker-Biocontainer-darkred?style=flat-square&maxAge=2678400)](https://quay.io/repository/biocontainers/bofasa?tab=info)

**B**acterial **O**rthology **F**inding **A**nd **S**yntenic **A**nalysis (**bofasa**)

<p align="center">
<img src="https://github.com/user-attachments/assets/11fa592f-656e-4578-93b8-17f8dc36f6c2" width="400">
</p>

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
- [Additional Scripts](#additional-scripts)
- [Documentation](#documentation)
- [Citation](#citation)
- [License](#license)

bofasa is specifically designed for investigating orthology between multiple-species of bacteria. It prioritizes high-quality orthology inference at the expense of throughput, designed to run on between 4 and 200 genomes.

For single species analyses, we recommend [Panaroo](https://github.com/gtonkinhill/panaroo) or [PPanGGOLiN](https://github.com/labgem/PPanGGOLiN), which for such cases, offers considerable advantages in terms of speed, accuracy, and scalability.

Comparable and good alternatives to bofasa include [PIRATE](https://github.com/SionBayliss/PIRATE) and [SCARAP](https://github.com/SWittouck/SCARAP). PIRATE performs similarly to bofasa when applied to genomes representative of different species from a single genus when using default parameters. SCARAP behaves more similar to non-bacteria specific multi-species orthology inference software such as OrthoFinder and SonicParanoid and is really fast and could be well-suited for more comprehensive and large-scale analysis. 

## Installation

### Bioconda

```bash
conda create -c conda-forge -c bioconda -p /path/to/bofasa_env/ bofasa
conda activate /path/to/bofasa_env/
bofasa setup

bofasa -h
```

### Pixi 

```bash
pixi init /path/to/bofasa_env/
pixi add -m /path/to/bofasa_env/ bofasa
pixi shell -m /path/to/bofasa_env/
bofasa setup

bofasa -h
```

### Docker (*via Biocondainters*)

*In the process of being developed.*

<!-- 
To get a Docker image from Quay.IO/Biocontainers, you can do something like the following, *note the platform designation might need to be adapted to your particular machine.* 

```bash
docker pull quay.io/biocontainers/bofasa:1.2.0--pyh106432d_0
```

Next, setup databases:

```bash
docker run -v /path/to/dbs/:/data/dbs/ -e BOFASA_DB_PATH=/data/dbs/ --platform linux/amd64 quay.io/biocontainers/bofasa:1.2.0--pyh106432d_0 bofasa setup
```

Here, `/path/to/dbs/` is the actual location on your computer where to store the databases and `/data/dbs/` is the location on the container it maps to. 

> [!IMPORTANT]
> When running `bofasa prep` and `bofasa run`, please issue the `-auto` flag to overcome interactive prompts. However, caution, this will lead to overwriting output directories!

> [!NOTE]
> You can also use Docker images via Singularity/Apptainer, can probably ask some AI agent how to do this. 
-->

### Test Installation:

```bash
# get test dataset and testing script from bofasa Github repo:
wget https://github.com/raufs/bofasa/raw/refs/heads/main/test_case.tar.gz
wget https://raw.githubusercontent.com/raufs/bofasa/refs/heads/main/run_tests.sh

bash run_tests.sh
```

> [!NOTE]
> `bofasa setup` should already have been run. Also, the above will work for conda, but dataset can be downloaded and adapted to test docker installation too!

## Quick Start

bofasa provides a unified command-line interface with two main subcommands: **`bofasa prep`** and **`bofasa run`**.

### 0. Setup databases (Pfam, geNomad, and ISFinder; *needs to be done only once!*)

```bash
bofasa setup
```

> [!NOTE]
> To aid reproducability between two runs, please ensure you are using the same version of the Pfam-A database!

### 1. Prepare Input Data

```bash
bofasa prep -i genome1.fasta genome2.fasta genome3.fasta genome4.fasta -o prep_output/
```

### 2. Run Analysis

```bash
bofasa run -i prep_output/ -o analysis_output/ -c 8
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
Processes input genomes (FASTA or GenBank files) for analysis. Can also take in [Prokka](https://github.com/tseemann/prokka) and [bakta](https://github.com/oschwengers/bakta) directories!

```bash
bofasa prep -i <genome1.fasta> <genome2.gbk> ... -o <output-dir> [OPTIONS]
```

**Required Arguments:**
- `-i, --input-genomes`: Input genomes (FASTA or GenBank files)
- `-a, --annotation-dirs`: Annotation directories (Prokka or Bakta output directories). Required if no input genomes are provided.
- `-o, --output-dir`: Output directory for prepared data

**Optional Arguments:**
- `-c, --threads`: Number of threads (default: 4)
- `-gcm, --gene-calling-method`: Gene calling method (pyrodigal/prodigal, default: pyrodigal)
- `-l, --locus-tag-length`: Length of locus tags (default: 3)
- `-m, --meta-mode`: Use meta-mode for gene calling
- `-ml, --min-length`: Minimum length of domain/inter-domain unit to consider (default: 20)
- `-rlt, --rename-locus-tags`: Rename locus tags in GenBank files and annotation directories (Prokka/Bakta)
- `-rg, --run-genomad`: Run genomad for phage/plasmid annotation
- `-emg, --extract-mge-genomes`: Extract mobile genetic element (phage/plasmid) genomes identified by genomad (requires -rg)
- `-sds, --skip-domain-splitting`: Skip domain-annotation based splitting of coding sequences
- `-mm, --max-memory`: Memory limit in GB (default: 32)
- `x, --ignore-upperbound-limit`: Ignore the upper bound limit for the number of genomes to process (default: 200)
- `-y, --auto`: Automatically set `y` for all interactive prompts

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
- `-dj, --dog-jaccard`: Jaccard index threshold (default: 0.5)
- `-fic, --fixation-index-cutoff`: Fixation index cutoff (default: 0.1)
- `-smb, --skip-merge-back`: Skip merge back assessment
- `-spr, --skip-phylo-refine`: Skip phylogenetic refinement
- `-rs, --rooting-seeds`: Rooting seeds (default: 1)
- `-us, --ultra-sens`: Use ultra-sensitivity mode
- `-mi, --mcl-inflation`: MCL inflation parameter (default: 1.2)
- `-ns, --near-scc-prop`: Near SCC proportion (default: 0.95)
- `-c, --threads`: Number of threads (default: 4)
- `-mrd, --max-recursion-depth`: Max recursion depth (default: 5000)
- `-mm, --max-memory`: Memory limit in GB (default: 32)
- `-y, --auto`: Automatically set `y` for all interactive prompts

## Additional Scripts

bofasa also provides several additional utility scripts:

- `extract_og_pairs.py` - Extract ortholog group pairs for method comparison
- `extract_scc_og_pairs.py` - Extract single-copy-core ortholog group pairs
- `compare_og_pairs.py` - Compare ortholog group pairs between methods
- `determine_og_freqs.py` - Determine ortholog group frequencies
- `print_og_itol_matrix.py` - Print ortholog group matrices for iTOL visualization
- `visualize_og_context.py` - Visualize the context of focal ortholog groups of interest. ***Is still experimental!***
  
**Note:** GenBank processing and Prodigal gene calling functionality is now integrated into the main `bofasa prep` command and no longer requires separate scripts.

## Documentation

For detailed information on bofasa, see the `ALGORITHM.md` document.

## Citation

> ***High-quality and automated inference of ortholog groups across multiple bacterial species using bofasa.***
> Rauf Salamzade, Eitan Yaffe, Aamuktha Kottapalli, Lindsay Kalan, David Relman, 2026.

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
