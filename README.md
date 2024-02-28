# bofasa
**B**acterial **O**rthology **F**inding **A**nd **S**yntenic **A**nalysis (**bofasa**)

_Developed by Rauf Salamzade, Aamuktha Kottapalli_

_Kalan lab @ University of Wisconsin-Madison; McMaster University_

bofasa is specifically designed for investigating orthology between multiple-species of bacteria and uses [OrthoFinder]() at its core to determine coarse ortholog groups. This is to leverage OrthoFinder's ability to normalize scores between proteins during orthology inference for both gene length discrepancies and genome-wide differences. Key adaptations are made however for appropriate application to bacterial genomes:

- 

For single species analyses, we recommend [panaroo](https://github.com/gtonkinhill/panaroo), which offers considerable advantages in terms of speed, accuracy, and scalability.

## Funding acknowledgment:

This project has been funded in whole or in part with Federal funds from the National Institute of Allergy and Infectious Diseases, National Institutes of Health.
