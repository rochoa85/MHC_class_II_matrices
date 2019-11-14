# MHC class II scoring-matrices

## Package to predict core region of peptides bound to MHC class II structure, model new peptides and run simulations

From publication *"Scoring-matrices from structural descriptors improve MHC class II motifs of bound peptides"*
Bioinformatics Journal (Oxford)
Authors: Rodrigo Ochoa, Roman Laskowski, Alessandro Laio, Janet Thornton, Pilar Cossio
Year: 2019

## Purpose

The goal of this script is to provide tools for applying scoring-matrices to predict core regions of peptides bound to MHC class structures. In addition, the user can model new peptides bound to crystal structures of MHC class II in complex with peptide templates to create their own matrices, and run sampling simulations using the backrub method from Rosetta

## Third-party tools required:

BioPython: https://biopython.org/wiki/Download
Rosetta Commons: https://www.rosettacommons.org/software/license-and-download

## How to run

One basic way to run the script is thorugh:

`model_peptides_MHC_complexes.py [-h] -l LIST_PEP [-m MODE] [-r ROSETTA]
                                       [-a ALLELE] [-o OUTPUT]`
                                       
The arguments are

```
optional arguments:
  -h, --help   show this help message and exit
  -l LIST_PEP  List with the peptides that want to be analyzed
  -m MODE      Choose a mode to run the script from thee options: 1) core, 2)
               model, 3) backrub
  -r ROSETTA   Version of Rosetta that will be implemented
  -a ALLELE    Allele of the MHC class II DRB1 that will be selected from four
               options: 1) 0101, 2) 0301, 3) 0401, 4) 1501
  -o OUTPUT    Name of the output file with the statistics results
 ```

## Example

`python3.5 model_peptides_MHC_complexes.py -l list_peptides.txt -m backrub`
