# The evolutionary dynamics of viral mimic-host protein interactions

This gitHub includes a collection of scripts and command-lines for running and reproducing the analysis workflow. The input data required for running the scripts described here can be found at our Zenodo dataset- https://zenodo.org/uploads/15880296

## Description

This repository contains code for data processing, statistical analysis, and visualization. It includes:

- **Command-lines for programs used in this study including flags and parameters used:
Guidance, 3SEQ, PAML, Selecton, BLAST, PhyML, AlphaFold, CSU, AlphaFold-Multimer, AIUpred, ScanNet


- **Detailed scripts for running complete analyses:
1) Structural homology analysis between human and viral proteins, matching viral mimicking - host mimicked proteins (FoldSeek and filtering of output results)
2) Positive selection analysis:
At the gene level, categorizing the studied proteins based on mimicry-related groups and testing for positive selection enrichment or depletion. 
At the residue level, categorizing residues by spatial classification (surface, interface, or core) and testing for positive selection enrichment or depletion.
3) Evolutionary rates analysis - Testing for higher\lower dN/dS in mimicry-related protein groups.
