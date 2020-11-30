# MHC class II scoring-matrices

## Package to predict core region of peptides bound to MHC class II structures, model new peptides and run simulations

* From publication *"Impact of structural interactions to predict the effect of single-point mutations in MHC class II peptide binders"*, 2021
* Authors: Rodrigo Ochoa, Roman A. Laskowski, Janet M. Thornton, Pilar Cossio

## Purpose

The goal of this script is to provide tools for applying scoring-matrices to predict core regions of peptides bound to MHC class structures. In addition, the user can model new peptides bound to crystal structures of MHC class II in complex with peptide templates to create their own matrices, and run sampling simulations using the backrub method from Rosetta

## Third-party tools required:

- BioPython: https://biopython.org/wiki/Download
- Rosetta Commons: https://www.rosettacommons.org/software/license-and-download

The BioPython module can be installed directly from package repositories. For the Rosetta functionalities the recommended is to follow the installation instructions, and take note of the Rosetta version that will be provided in the script.

## How to run

The basic command line to run the script is:

`model_peptides_MHC_complexes.py [-h] -l LIST_PEP [-m MODE] [-r ROSETTA]
                                       [-a ALLELE] [-o OUTPUT]`
                                       
where the arguments are:

```
optional arguments:
  -h, --help   show this help message and exit
  -l LIST_PEP  List with the peptides that want to be analyzed (required)
  -m MODE      Choose a mode to run the script from thee options: 1) core, 2)
               model, 3) backrub
  -r ROSETTA   Version of Rosetta that will be implemented
  -a ALLELE    Allele of the MHC class II DRB1 that will be selected from four
               options: 1) 0101 (default), 2) 0301, 3) 0401, 4) 1501
  -o OUTPUT    Name of the output file with the statistics results
 ```
The only required argument is the file containing the peptides. For the other parameters the script has default values. **However, please be aware of changing the Rosetta version through the flag or directly in the script by default.**

## Examples

### Predict core regions

In case the user want to predict the core region of several peptides, the way to run the script is:

`python3.5 model_peptides_MHC_complexes.py -l list_peptides.txt -m core`

After that, the results can be seen in the screen as (the results are also saved in an output file):
```
Original peptide sequence: IPTAFSIGKTYKPEE
Predicted core with structural descriptors:  FSIGKTYKP
Predicted core with sequence-based motifs:  SIGKTYKPE
```

### Model peptides bound to a particular allele

To model peptides of interest in an allele with crystal structure available, the user can select from four alleles that have available crystals in the local folder provided. Based on the selection, the script reads the file `list_MHC_crystals.txt` and select the corresponding structure. One way to run the script is as follows:

`python3.5 model_peptides_MHC_complexes.py -l list_peptides.txt -m model -a 0101`

The file `<sequence>_modelled.pdb` will be generated. In addition, a report with the statistics of the modelling will be created. By default the name of the file is `stats_peptide_models.txt`- However, the name can be changed with the `-o` flag. An example of the report is the following:

```
Crystal	Allele	Peptide_template	Peptide_model	Core_template	Core_model	Alignment_score	Identities	Identical_pockets
1t5x	DRB1*01:01	AAYSDQATPLLLSPR	IPTAFSIGKTYKPEE	FSIGKTYKP	YSDQATPLL	0	2	1
The model was completed
```
In addition to the prediction of the cores, the peptides are aligned to evaluate how similar they are, and identities are mapped in the full sequence and in the core region. This is useful to understand how different is the sequence content of the peptides from the available peptide template.

### Run additional sampling using Rosetta Backrub

After modelling the peptides, it is possible to run a simulation of the system using the backrub method from Rosetta. For that purpose, just run the script with the following options:

`python3.5 model_peptides_MHC_complexes.py -l list_peptides.txt -m backrub -a 0101`

In that case, the model with the additional simulation will be run in the local computer. The output of the trajectory can be used to calculate any structural descriptor of interest using available packages to study interactions, or calculate other variables such as secondary structures or accessible surface areas, among others.

## Support

In case the protocol is useful for other research projects and require some advice, please contact us to the email: rodrigo.ochoa@udea.edu.co
