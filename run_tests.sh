#!/usr/bin/env bash

# Step 0: Uncompress test_case.tar.gz and cd into it.
rm -rf test_case/
tar -zxvf test_case.tar.gz 
cd test_case/

# run bofasa on test set of four Staphylococcus genomes
# Assumes bofasa setup is already run!!!
bofasa prep -i Genomes/*.fna -o Bofasa_Prep_Results/ -c 4
bofasa run -i Bofasa_Prep_Results/ -o Bofasa_Run_Results/ -c 4