#!/usr/bin/python

"""
Package to model peptides bound to MHC class II complexes and run backrub simulations

From publication "Scoring-matrices from structural descriptors improve MHC class II motifs of bound peptides"
Journal: Oxford Bioinformatics 
Authors: Rodrigo Ochoa, Roman Laskowski, Alessandro Laio, Janet Thornton, Pilar Cossio
Year: 2019

Third-party tools required:

BioPython: https://biopython.org/wiki/Download
Rosetta Commons: https://www.rosettacommons.org/software/license-and-download
"""

########################################################################################
# Authorship
########################################################################################

__author__ = "Rodrigo Ochoa"
__credits__ = ["Rodrigo Ochoa", "Roman Laskowski", "Alessandro Laio","Janet Thornton","Pilar Cossio"]
__license__ = "GPL"
__version__ = "1.0"
__email__ = "rodrigo.ochoa@udea.edu.co"

########################################################################################
# Modules to import
########################################################################################

import os
import argparse
import subprocess
from random import shuffle

# Third-party modules
from Bio.PDB import *
from Bio import pairwise2

########################################################################################
# Functions
########################################################################################

def replace_amino_acid(pep_pdb,pep_chain,old_aa,new_aa,pep_position):
    """
    Left only backbone atoms for the replacement of an amino acid in a position of the pdb chain
        
    Arguments:
    pep_pdb -- Biopython object with the pdb structure of interest
    pep_chain -- Chain containing the sequence where the replacement will be done
    old_aa -- Amino acid that will be replaced
    new_aa -- Amino acid that will be mutated
    pep_position -- Position in the chain where the mutation will be done
    
    Output:
    A PDB structure called "pre-mutated.pdb" with only the backbone atoms
    """
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}

    # Read the PDB file
    residues=pep_pdb.get_residues()
    chain=pep_pdb[0][pep_chain]
    
    # Report the mutation made
    message="The residue {} in chain {} and position {} will be changed by {}".format(old_aa,pep_chain,pep_position,new_aa)
    print(message)
    
    # Rename the residue
    chain[pep_position].resname=aminoacids[new_aa]
    
    # Delete the other atoms leaving only the atoms of the backbone
    ids=[]
    for a in chain[pep_position]:
        atomId=a.id
        if atomId not in ("N","CA","O","C"): ids.append(atomId)
    for i in ids: chain[pep_position].detach_child(i)
        
    # Saving the new structure
    io = PDBIO()
    io.set_structure(pep_pdb)
    io.save("pre-mutated.pdb")

##################################################################################

def generatePairs(pep1,pep2):
    """
    Definition of the pair of amino acids that will be replaced with each other between 2 peptide sequences
    
    Arguments:
    pep1 -- First peptide sequence
    pep2 -- Second peptide sequence
    
    Return:
    mutation_list -- List containing tuples with the old and new amino acids, as well as their positions
    """
    
    mutation_list=[]
    
    # Iteration over the list of amino acids of the petide sequences
    for i,aa in enumerate(pep1):
        mutation_list.append((aa,pep2[i],i+1))
    
    # Return the list of tuples
    return mutation_list

##################################################################################

def generateMutations(pep1,pep2):
    """
    Definition of the pair of amino acids that will be replaced with each other between 2 peptide sequences, but taking into account
    insertions or avoiding the replacement if the amino acid is the same
    
    Arguments:
    pep1 -- First peptide sequence
    pep2 -- Second peptide sequence
    
    Return:
    mutation_list -- List containing tuples with the old and new amino acids, as well as their positions
    """
    mutation_list=[]
    
    # Check if the position is open to be inserted an amino acid
    if pep1[0]=="X":
        counter=0
        for i,aa in enumerate(pep1):
            if aa=="X": pass
            else:
                counter+=1
                # Check that both amino acids are different
                if aa!=pep2[i]:
                    mutation_list.append((aa,pep2[i],counter))
    else:
        # Do the normal process, but verifyng that both amino acids are different
        for i,aa in enumerate(pep1):
            if aa!=pep2[i]:
                mutation_list.append((aa,pep2[i],i+1))
    
    # Return the list of tuples
    return mutation_list

##################################################################################

def mutateRosetta(pep_chain,new_aa,pep_position,rosetta_path):
    """
    Prediction of the mutated side chain using Rosetta functionalities
    
    Arguments:
    pep_chain -- Chain containing the sequence where the replacement will be done
    new_aa -- Amino acid that will be mutated
    pep_position -- Position in the chain where the mutation will be done
    rosetta_path -- Path of the Rosetta installation in the local computer
    
    Output:
    A PDB structure called "post-mutated.pdb" with the predicted side chain atoms
    """
    
    # Generate Rosetta config file with instructions to perform the mutation
    rosetta_config_file=open("resfile.config","w")
    rosetta_config_file.write("NATRO\n")
    rosetta_config_file.write("start\n")
    rosetta_config_file.write("\t{} {} PIKAA {} EX 1 EX 2 EX 3 EX 4 EX_CUTOFF 0\n".format(pep_position,pep_chain,new_aa))
    rosetta_config_file.close()
    
    # Run in a BASH environment the Fixbb tool from Rosetta
    os.system("{}/main/source/bin/fixbb.default.linuxgccrelease -in:file:s pre-mutated.pdb -resfile resfile.config".format(rosetta_path))
    # Open the structure with Biopython and save a clean pdb file
    parser = PDBParser()
    structure = parser.get_structure('REF',"pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-pre-mutated.pdb")
    
    # Run in a BASH environment the Relax protocol from Rosetta with backbone fixed
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s post-pre-mutated.pdb -relax:thorough -relax:bb_move false".format(rosetta=rosetta_path))
    # Open the structure with Biopython and save a clean pdb file
    parser = PDBParser()
    structure = parser.get_structure('REF',"post-pre-mutated_0001.pdb")
    io = PDBIO()
    io.set_structure(structure)
    io.save("post-mutated.pdb")
    
################################################################

def get_blueprint_add(mutList,numberPrint):
    """
    Generation of a blueprint required for the Rosetta Remodel package to add a single amino acids in any flanking region
    
    Arguments:
    mutList -- List of the mutations that will be performed, to prioritize the insertion of new amino acids
    numberPrint -- Number assigned to the blueprint file (in this case the value will be 1)
    
    Output:
    A blueprint file with the instructions to perform the amino acid insertion
    """
    
    # Opening the blueprint file
    print("Processing blueprint {} ...".format(numberPrint))
    reportFile=open("blueprint{}.txt".format(numberPrint),"w")
    
    # Iterate over the list of mutaitons
    for mut in mutList:
        old_aa=mut[0]
        new_aa=mut[1]
        position=mut[2]
        
        # Check if the amino acid should be inserted, if should relax the amino acid next to it, or if the amino acid will remain fixed
        if new_aa=="X":
            reportFile.write("{} {} .\n".format(position,old_aa))
        elif position==len(mutList)-1:
            reportFile.write("{} {} L PIKAA {}\n".format(len(mutList)-1,new_aa,new_aa))
        elif position==len(mutList):
            reportFile.write("0 X L PIKAA %s EX 1 EX 2 EX 3 EX 4 EX_CUTOFF 0\n" %new_aa)
        else:
            reportFile.write("{} {} .\n".format(position,new_aa))
    
    # Close the file with the Rosetta Remodel instructions
    reportFile.close()

################################################################

def get_blueprint_add2(mutList,numberPrint):
    """
    Generation of a blueprint required for the Rosetta Remodel package to add a second amino acid after adding one in any flanking region
    
    Arguments:
    mutList -- List of the mutations that will be performed, to prioritize the insertion of new amino acids
    numberPrint -- Number assigned to the blueprint file (in this case the value will be 1)
    
    Output:
    A blueprint file with the instructions to perform the amino acid insertion
    """
    
    # Opening the blueprint file
    print("Processing blueprint {} ...".format(numberPrint))
    reportFile=open("blueprint{}.txt".format(numberPrint),"w")
    
    # Iterate over the list of mutaitons
    for mut in mutList:
        old_aa=mut[0]
        new_aa=mut[1]
        position=mut[2]
        
        # Check if the amino acid should be inserted, if should relax the amino acid next to it, or if the amino acid will remain fixed
        if new_aa=="X":
            reportFile.write("{} {} .\n".format(position-int(mutList[0][2]),old_aa))
        elif position==int(mutList[0][2])+1:
            reportFile.write("{} {} L PIKAA {}\n".format(1,new_aa,new_aa))
        elif position==int(mutList[0][2]):
            reportFile.write("0 X L PIKAA %s EX 1 EX 2 EX 3 EX 4 EX_CUTOFF 0\n" %new_aa)
        else:
            reportFile.write("{} {} .\n".format(position-int(mutList[0][2]),new_aa))
    
    # Close the file with the Rosetta Remodel instructions
    reportFile.close()
    
#################################################################

def get_chains(pdb_file):
    """
    Auxiliary function to extract all the chains from a structure in separate PDB files
    
    Arguments:
    pdb_file -- Crystal structure with the chains of interest
    
    Output:
    The PDB files with the corresponding chain letter in the name of the file
    """
    
    # Read the crystal structure as a Biopython object
    parser = PDBParser()
    structure = parser.get_structure('PEP',pdb_file)
    
    # Iterate over the chains
    for chain in structure[0]:
        ppb = CaPPBuilder()
        chLetter=chain.get_id()
        name=pdb_file.split(".")[0]
        
        # Store the sequence in case it is required
        for pp in ppb.build_peptides(chain):
            seq=pp.get_sequence()
        
        # Save the PDB file
        io = PDBIO()
        io.set_structure(chain)
        io.save('{}_{}.pdb'.format(name,chLetter))
        
#################################################################

def remodel_single_aa(pep_to_model,rosetta_path,old_aa,position):
    """
    Run the Rosetta Remodel protocol for inserting the first and single amino acid in any flanking region
    NOTE: This function depends on the availability of a BASH environment - Tested on Ubuntu 16.04
    
    Arguments:
    pep_to_model -- Full sequence of the peptide that is going to be modelled
    rosetta_path -- Path of the Rosetta installation in the local computer
    old_aa -- Amino acid that will be replaced
    position -- Position in the sequence of the amino acid that will be replaced
    
    Output:
    A PDB structure of the peptide modelled in complex to the MHC class II receptor allele 
    """
    
    # List of amino acids
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
    
    # Split the chains
    get_chains("{}_pre_mutation.pdb".format(pep_to_model))
    
    # Save the chains in the files required to run the modelling
    os.system("cat {pep}_pre_mutation_A.pdb {pep}_pre_mutation_B.pdb | grep -v END > target.pdb".format(pep=pep_to_model))
    os.system("mv {}_pre_mutation_C.pdb only_chain.pdb".format(pep_to_model))
    os.system("rm {pep}_pre_mutation_A.pdb {pep}_pre_mutation_B.pdb".format(pep=pep_to_model))
    
    # Run the Remodel protocol and save the result in a separate file
    os.system("{rosetta}/main/source/bin/remodel.default.linuxgccrelease -database {rosetta}/main/database -in:file:s only_chain.pdb -remodel:blueprint blueprint1.txt -run:chain C -overwrite".format(rosetta=rosetta_path))
    os.system("csplit only_chain_0001.pdb /All/")
    os.system("mv xx00 peptide_pre.pdb")
    os.system("sed -i 's# A  # C  #g' peptide_pre.pdb")
    os.system("rm 1.pdb 2.pdb 3.pdb 4.pdb 5.pdb xx01 only_chain*")
    
    # Deletion of the amino acid from the template that was not changed
    if len(str(position))==1:
        os.system("grep -v '{} C   {}' peptide_pre.pdb > peptide_pre2.pdb".format(aminoacids[old_aa],position))
    if len(str(position))==2:
        os.system("grep -v '{} C  {}' peptide_pre.pdb > peptide_pre2.pdb".format(aminoacids[old_aa],position))
    
    # Generate Rosetta config file
    rosetta_config_file=open("resfile.txt","w")
    rosetta_config_file.write("NATRO\n")
    rosetta_config_file.write("start\n")
    rosetta_config_file.close()
    
    # Renumber of the PDB amino acids
    os.system("{}/main/source/bin/fixbb.default.linuxgccrelease -in:file:s peptide_pre2.pdb -resfile resfile.txt -renumber_pdb -per_chain_renumbering".format(rosetta_path))
    os.system("sed -i 's# A  # C  #g' peptide_pre2_0001.pdb")
    os.system("csplit peptide_pre2_0001.pdb /All/")
    os.system("mv xx00 peptide.pdb; rm peptide_pre* xx01")
    
    # Final relaxation of the system
    os.system("cat target.pdb peptide.pdb > prev.pdb")
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s prev.pdb -relax:fast -relax:bb_move false".format(rosetta=rosetta_path))
    os.system("rm target.pdb peptide.pdb")
    os.system("csplit prev_0001.pdb /All/")
    os.system("mv xx00 {}_modelled.pdb".format(pep_to_model))
    os.system("rm xx01 prev_0001.pdb score.sc prev.pdb {}_pre_mutation.pdb".format(pep_to_model))

#################################################################

def remodel_double_aa(pep_to_model,rosetta_path,old_aa_1,position_1,old_aa_2,position_2):
    """
    Run the Rosetta Remodel protocol for inserting two amino acids in any flanking region
    NOTE: This function depends on the availability of a BASH environment - Tested on Ubuntu 16.04
    
    Arguments:
    pep_to_model -- Full sequence of the peptide that is going to be modelled
    rosetta_path -- Path of the Rosetta installation in the local computer
    old_aa_1 -- First amino acid that will be replaced
    position_1 -- Position in the sequence of the first amino acid that will be replaced
    old_aa_2 -- Second amino acid that will be replaced
    position_2 -- Position in the sequence of the second amino acid that will be replaced
    
    Output:
    A PDB structure of the peptide modelled in complex to the MHC class II receptor allele 
    """
    
    aminoacids={"A":"ALA","D":"ASP","E":"GLU","F":"PHE","H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","G":"GLY",
                "N":"ASN","P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL","W":"TRP","Y":"TYR","C":"CYS"}
    
    # Split the chains and perform the mutation
    get_chains("{}_pre_mutation.pdb".format(pep_to_model))
    
    # Save the chains in the files required to run the modelling
    os.system("cat {pep}_pre_mutation_A.pdb {pep}_pre_mutation_B.pdb | grep -v END > target.pdb".format(pep=pep_to_model))
    os.system("mv {}_pre_mutation_C.pdb only_chain.pdb".format(pep_to_model))
    os.system("rm {pep}_pre_mutation_A.pdb {pep}_pre_mutation_B.pdb".format(pep=pep_to_model))
    
    # Run the Remodel protocol for the first time and save the result in a separate file
    os.system("{rosetta}/main/source/bin/remodel.default.linuxgccrelease -database {rosetta}/main/database -in:file:s only_chain.pdb -remodel:blueprint blueprint1.txt -run:chain C -overwrite".format(rosetta=rosetta_path))
    os.system("sed -i 's# A  # C  #g' only_chain_0001.pdb")
    os.system("rm 1.pdb 2.pdb 3.pdb 4.pdb 5.pdb")
    
    # Run the Remodel protocol for the second time and save the result in a separate file
    os.system("{rosetta}/main/source/bin/remodel.default.linuxgccrelease -database {rosetta}/main/database -in:file:s only_chain_0001.pdb -remodel:blueprint blueprint2.txt -run:chain C -overwrite".format(rosetta=rosetta_path))
    os.system("csplit only_chain_0001_0001.pdb /All/")
    os.system("mv xx00 peptide_pre.pdb")
    os.system("sed -i 's# A  # C  #g' peptide_pre.pdb")
    os.system("rm 1.pdb 2.pdb 3.pdb 4.pdb 5.pdb xx01 only_chain*")
    
    # Deletion of the amino acid from the template that was not changed
    if len(str(position_1))==1:
        os.system("grep -v '{} C   {}' peptide_pre.pdb | grep -v '{} C   {}' > peptide_pre2.pdb".format(aminoacids[old_aa_1],position_1,aminoacids[old_aa_2],position_2))
    if len(str(position_1))==2:
        os.system("grep -v '{} C  {}' peptide_pre.pdb | grep -v '{} C  {}' > peptide_pre2.pdb".format(aminoacids[old_aa_1],position_1,aminoacids[old_aa_2],position_2))
    
    # Generate Rosetta config file
    rosetta_config_file=open("resfile.txt","w")
    rosetta_config_file.write("NATRO\n")
    rosetta_config_file.write("start\n")
    rosetta_config_file.close()
    
    # Renumber of the PDB amino acids
    os.system("{}/main/source/bin/fixbb.default.linuxgccrelease -in:file:s peptide_pre2.pdb -resfile resfile.txt -renumber_pdb -per_chain_renumbering".format(rosetta_path))
    os.system("sed -i 's# A  # C  #g' peptide_pre2_0001.pdb")
    os.system("csplit peptide_pre2_0001.pdb /All/")
    os.system("mv xx00 peptide.pdb; rm peptide_pre* xx01")
    
    # Final relaxation of the system
    os.system("cat target.pdb peptide.pdb > prev.pdb")
    os.system("{rosetta}/main/source/bin/relax.default.linuxgccrelease -database {rosetta}/main/database -in:file:s prev.pdb -relax:fast -relax:bb_move false".format(rosetta=rosetta_path))
    os.system("rm target.pdb peptide.pdb")
    os.system("csplit prev_0001.pdb /All/")
    os.system("mv xx00 {}_modelled.pdb".format(pep_to_model))
    os.system("rm xx01 prev_0001.pdb score.sc prev.pdb {}_pre_mutation.pdb".format(pep_to_model))

#################################################################

def check_motif_seq(sequence,allele):
    """
    Prediction of the potential core regions based on available sequence-based motifs
    source: http://www.cbs.dtu.dk/biotools/MHCMotifViewer/Home.html
        
    Arguments:
    sequence -- Peptide which core will be predicted
    allele -- MHC class II allele that will be used to obtain the scores
        
    Return:
    core -- The fragment with the highest probability of being the core
    """
    
    # Read the matrix of the chosen allele
    if allele=="1501":
        matrix=[x.strip() for x in open("matrices/motif_seq_DRB1_1501.txt")]
    if allele=="0101":
        matrix=[x.strip() for x in open("matrices/motif_seq_DRB1_0101.txt")]
    if allele=="0301":
        matrix=[x.strip() for x in open("matrices/motif_seq_DRB1_0301.txt")]
    if allele=="0401":
        matrix=[x.strip() for x in open("matrices/motif_seq_DRB1_0401.txt")]

    # Order of the amino acids in the original matrix
    listAA=["A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

    #Store the matrix data in a dictionary per amino acid and position
    matrixScores={}
    for i,m in enumerate(matrix):
        data=m.split()
        for j,d in enumerate(data):
            matrixScores[listAA[j]+"-"+str(i+1)]=d
    
    # Read the sequence by fragments of 9 AA and calculate the core values
    coreScores={}
    for i in range(9,len(sequence)+1):
        fragment=sequence[i-9:i]
        score_core=0
        for pos,aa in enumerate(fragment):
            score_core+=float(matrixScores[aa+"-"+str(pos+1)])
        coreScores[fragment]=score_core
     
    # Choose the fragment with the HIGHEST score according to the motifs definition
    core=min(coreScores, key=coreScores.get)
    
    # Return the fragment predicted as the core 
    return core

#################################################################

def check_descriptor_struct(sequence,descriptor="contacts"):
    """
    Prediction of the potential core regions based on the calculated structural descriptors
        
    Arguments:
    sequence -- Peptide which core will be predicted
    descriptor -- There are two options: "contacts" (default) or "hbonds_mc"
        
    Return:
    core -- The fragment with the highest probability of being the core
    """
    
    # Read file with the descriptor chosen
    matrix=[x.strip() for x in open("matrices/descriptor_{descriptor}_DRB1.txt".format(descriptor=descriptor))]
    
    #Store the matrix data in a dictionary per amino acid and position
    matrixScores={}
    for i,m in enumerate(matrix):
        if i!=0:
            data=m.split()
            for j,d in enumerate(data):
                if j==0: aa=d
                else: matrixScores[aa+"-"+str(j)]=d
    
    # Read the sequence by fragments of 9 AA and calculate the core values
    coreScores={}
    for i in range(9,len(sequence)+1):
        fragment=sequence[i-9:i]
        score_core=0
        for pos,aa in enumerate(fragment):
            score_core+=float(matrixScores[aa+"-"+str(pos+1)])
        coreScores[fragment]=score_core
     
    # Choose the fragment with the LOWEST score according to the descriptor definition
    core=min(coreScores, key=coreScores.get)
    
    # Return the fragment predicted as the core 
    return core

#################################################################

def obtain_peptides_target_template(core_final,crystals,statsFile):
    """
    Choose the peptide template from the available crystal files and save the report
        
    Arguments:
    core_final -- Core sequence predicted for the peptide of interest
    crystals -- List with the information of the available crystals
    statsFile -- Object with the file that will print the report
        
    Return:
    pep_template -- Final sequence of the template peptide,
    pep_model -- Final sequence of the peptide to model
    crystal -- Crystal used to obtain the structure
    """
    
    # Read across all the available crystals
    for rows in crystals:
        info=rows.split()
        allele_to_compare=info[1]
        
        # Check if the allele to compare is the same of the crystal
        if allele==allele_to_compare:
            # Alignment of the core region and calculation of the score
            crystal=info[0]
            peptide_template=info[5] # With X is represented missing or new amino acids
            alignment = pairwise2.align.globalxx(peptide_template,core_final)
            align_flag=0
            for a in alignment:
                if a[4]==len(peptide_template):
                    align_flag=1
                    score=a[2]
            if align_flag==0: score=0
                                
            # Check the conservation of pockets
            count_identities=0
            count_pocket=0
            for j,c_pos in enumerate(core_final):
                if c_pos==peptide_template[j]: count_identities+=1
                # Check for key amino acids
                if j==0:
                    if c_pos==peptide_template[j]: count_pocket+=1
                if j==3:
                    if c_pos==peptide_template[j]: count_pocket+=1
                if j==5:
                    if c_pos==peptide_template[j]: count_pocket+=1
                if j==8:
                    if c_pos==peptide_template[j]: count_pocket+=1
                    
            # Check how many AA were located in the flanking regions
            full_pep_template=info[4]
            
            # Get the flanking regions for the template
            for z in range(9,len(full_pep_template)+1):
                fragment=full_pep_template[z-9:z]
                if fragment==peptide_template:
                    nTer_templ=full_pep_template[0:z-9]
                    cTer_templ=full_pep_template[z:]
                    
            # Get the flanking regions for the model
            for z in range(9,len(pep_to_model)+1):
                fragment=pep_to_model[z-9:z]
                if fragment==core_final:
                    nTer_model=pep_to_model[0:z-9]
                    cTer_model=pep_to_model[z:]
            
            # Check the differences between the 2 peptides
            pep_template=""
            pep_model=""
            
            # If the peptide sizes are the same (based on the N-terminal region)
            if len(nTer_model)==len(nTer_templ):
                pep_template=full_pep_template
                pep_model=pep_to_model
                
            # if the N-terminal region of the peptide is longer
            elif len(nTer_model)>len(nTer_templ):
                diff=abs(len(nTer_model)-len(nTer_templ))
                pep_template="X"*diff+nTer_templ+peptide_template+cTer_templ
                pep_model=nTer_model+core_final+cTer_model+"X"*diff
            
            # if the N-terminal region of the peptide is shorter
            elif len(nTer_templ)>len(nTer_model):
                diff=abs(len(nTer_templ)-len(nTer_model))
                pep_template=nTer_templ+peptide_template+cTer_templ+"X"*diff
                pep_model="X"*diff+nTer_model+core_final+cTer_model
            
            #print info[0],info[1],full_pep_template,core_final,peptide_template,score,count_identities,count_pocket
            statsFile.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(info[0],info[1],full_pep_template,pep_to_model,core_final,peptide_template,score,count_identities,count_pocket))
            
            # Return the final sequences of the template peptide, the peptide to model and the crystal used to obtain the structure
            return pep_template,pep_model,crystal

#################################################################

def modelling(pep1_final,pep2_final,pep_to_model,crystal,rosetta_version):
    """
    Modelling and insertion of amino acids of the peptide to model
        
    Arguments:
    pep1_final -- Final sequence of the peptide used as template
    pep2_final -- Final sequence of the peptide to model
    pep_to_model -- Original sequence of the peptide to model
    crystal -- Crystal structure used to obtain the peptide template
    rosetta_version -- Version of Rosetta to implement in the protocol
        
    Output:
    PDB file of the peptide modelled based on the requirements
    """
    
    print("Starting modelling of peptide {} ...".format(pep_to_model))
    # Define the peptide chain letter
    pepChain="C"
    
    # Get the rosetta path
    bash = "locate -b {} | head -n1".format(rosetta_version)
    rosetta_path = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    
    # Generate the pairs of amino acids that will be mutated
    mutList=generateMutations(pep1_final,pep2_final)
    
    # Print and shuffle the list of mutations
    print(mutList)
    shuffle(mutList)
    
    # Iterate over the list of required mutations
    counter_step=0
    for j,mut in enumerate(mutList):
        old_aa=mut[0]
        new_aa=mut[1]
        position=mut[2]
        
        # Mutate positions with corresponding amino acids
        if old_aa != "X" and new_aa != "X":
            if counter_step==0:
                # Copy the crystal chosen to run the modelling in the first step and load a Biopython object
                os.system("cp crystals/{}.pdb .".format(crystal))
                parser = PDBParser()
                reference = parser.get_structure('REF',"{}.pdb".format(crystal))
                counter_step+=1
            else:
                # Load the previous mutated structure in a Biopython object
                parser = PDBParser()
                reference = parser.get_structure('REF',"post-mutated.pdb")
            
            # Replace the amino acid
            replace_amino_acid(reference,pepChain,old_aa,new_aa,position)
            mutateRosetta(pepChain,new_aa,position,rosetta_path)
            os.system("rm pre-mutated.pdb pre-mutated_0001.pdb post-pre-mutated.pdb post-pre-mutated_0001.pdb resfile.config score.sc")
    
    # Remove files not required for the insertions
    os.system("mv post-mutated.pdb {}_pre_mutation.pdb".format(pep_to_model))
    os.system("rm {}.pdb".format(crystal))
    
    # Check if the structure do not need of insertions
    if aa_to_model==0:
        os.system("cp {pep}_pre_mutation.pdb {pep}_modelled.pdb".format(pep=pep_to_model))
    
    # If it is required to add only one amino acid
    if aa_to_model==1:
        # Generate all the pairs of amino acids
        pairList=generatePairs(pep1_final,pep2_final)
        
        # Check if the amino acid is an insertion
        if pairList[0][1]=="X":
            # Generate the corresponding blueprint
            get_blueprint_add(pairList,1)
        else:
            # Generate the corresponding blueprint
            get_blueprint_add2(pairList,1)
           
        for pair in pairList:
            old_aa=pair[0]
            new_aa=pair[1]
            position=pair[2]
            
            # Remodel if this is the amino acid that need to be deleted
            if new_aa=="X":
                remodel_single_aa(pep_to_model,rosetta_path,old_aa,position)
                os.system("rm blueprint*")
                
    # If it is required to add two amino acids
    elif aa_to_model==2:
        # Generate all the pairs of amino acids
        pairList=generatePairs(pep1_final,pep2_final)
        
        # Check if the amino acid is an insertion
        if pairList[0][1]=="X":
            # Generating the corresponding blueprints
            get_blueprint_add(pairList[:-1],1)
            get_blueprint_add(pairList,2)
        else:
            # Generating the corresponding blueprints
            get_blueprint_add2(pairList[1:],1)
            get_blueprint_add2(pairList,2)

        # List containing the amino acids to include
        aa_list=[]
        position_list=[]
        for pair in pairList:
           old_aa=pair[0]
           new_aa=pair[1]
           position=pair[2]
            
           if new_aa=="X":
               aa_list.append(old_aa)
               position_list.append(position)
        
        # Remodel if these are the amino acids that need to be deleted
        remodel_double_aa(pep_to_model,rosetta_path,aa_list[0],position_list[0],aa_list[1],position_list[1])
        os.system("rm blueprint*")
        
#################################################################

def run_backrub(pep_to_model,rosetta_version):
    """
    Function to run Backrub method from Rosetta
    
    Arguments:
    pep_to_model -- Sequence of the peptide that was modelled and requires sampling
    rosetta_version -- Version of Rosetta that will be used for runnning backrub
    
    Output:
    
    """
    
    # Get the rosetta path
    bash = "locate -b {} | head -n1".format(rosetta_version)
    rosetta_path = subprocess.check_output(['bash','-c', bash]).strip().decode("utf-8")
    
    # Run backrub with default parameters
    os.system("{rosetta}/main/source/bin/backrub.default.linuxgccrelease -database {rosetta}/main/database \
              -s {pep}_modelled.pdb -ignore_unrecognized_res -ex1 -ex2 -extrachi_cutoff 0 -backrub:ntrials 2000 -mc_kt 1.2 \
              -ignore_zero_occupancy=false -trajectory=true -initial_pack".format(rosetta=rosetta_path,pep=pep_to_model))
    
    # Rename the files
    os.system("mv {pep}_modelled_0001_last.pdb {pep}_modelled_backrub_last.pdb".format(pep=pep_to_model))
    os.system("mv {pep}_modelled_0001_low.pdb {pep}_modelled_backrub_low.pdb".format(pep=pep_to_model))
    os.system("mv {pep}_modelled_0001.pdb {pep}_modelled_backrub.pdb".format(pep=pep_to_model))
    os.system("mv {pep}_modelled_0001_traj.pdb {pep}_modelled_backrub_traj.pdb".format(pep=pep_to_model))
    os.system("rm score.sc")

########################################################################################
########################################################################################
########################################################################################
# Main execution
########################################################################################
########################################################################################
########################################################################################
if __name__ == '__main__':
    
    # Script arguments
    parser = argparse.ArgumentParser(description='Prediction of core regions of peptides bound to MHC class II and modelling of peptides based on crystal template')
    parser.add_argument('-l', dest='list_pep', action='store',required=True,
                        help='List with the peptides that want to be analyzed')
    parser.add_argument('-m', dest='mode', action='store', default="core",
                        help='Choose a mode to run the script from thee options: 1) core, 2) model, 3) backrub')
    parser.add_argument('-r', dest='rosetta', action='store', default="rosetta_src_2016.32.58837_bundle",
                        help='Version of Rosetta that will be implemented')
    parser.add_argument('-a', dest='allele', action='store', default="0101",
                        help='Allele of the MHC class II DRB1 that will be selected from four options: 1) 0101, 2) 0301, 3) 0401, 4) 1501')
    parser.add_argument('-o', dest='output', action='store', default="stats_peptide_models.txt",
                        help='Name of the output file with the statistics results')
    args = parser.parse_args()
    
    ####################################################################################
    # Assignment of parameters
    
    # Rosetta version installed in the computer
    rosetta_version=args.rosetta
    
    # Allele that will be used to look for a template
    possible_alleles=["0101","0301","0401","1501"]
    if args.allele in possible_alleles:
        allele_short=args.allele
        if allele_short=="0101": allele="DRB1*01:01"
        if allele_short=="0301": allele="DRB1*03:01"
        if allele_short=="0401": allele="DRB1*04:01"
        if allele_short=="1501": allele="DRB1*15:01"
    else:
        print("The required allele is invalid. Please check the help option [-h]")
        quit()
    
    # List of peptides to model
    list_peptides=[x.strip() for x in open(args.list_pep)]
    
    # Output files
    statsFile=open(args.output,"w")
    
    # Check the mode
    possible_modes=["core","model","backrub"]
    if args.mode not in possible_modes:
        print("The required mode is invalid. Please check the help option [-h]")
        quit()
        
    ####################################################################################
    
    # Information crystals
    crystals=[x.strip() for x in open("list_MHC_crystals.txt")]
    
    # Iterate over the list of peptides
    for pep_to_model in list_peptides:
        
        # 1. Predict the cores using one of the methodologies    
        print("Original peptid sequence: {}".format(pep_to_model))
        
        core_final=check_descriptor_struct(pep_to_model,"contacts")
        print("Predicted core with structural descriptors: ",core_final)
        
        core_motif=check_motif_seq(pep_to_model,allele_short)
        print("Predicted core with sequence-based motifs: ",core_motif)
        print()
        
        if args.mode=="model" or args.mode=="backrub":
            # 2. Identification of the template and model peptides
            statsFile.write("Crystal\tAllele\tPeptide_template\tPeptide_model\tCore_template\tCore_model\tAlignment_score\tIdentities\tIdentical_pockets\n")
            pep_template,pep_model,crystal=obtain_peptides_target_template(core_final,crystals,statsFile)
            print(pep_template,pep_model)
                        
            # 3. Check the criterium of number of amino acids required to model to generate the peptide
            aa_to_model=abs(len(pep_template)-len(pep_to_model))
                        
            if aa_to_model<=2:
                modelling(pep_template,pep_model,pep_to_model,crystal,rosetta_version)
                statsFile.write("The model was completed\n")
                # Run the Backrub simulation in case was required
                if args.mode=="backrub":
                    run_backrub(pep_to_model,rosetta_version)
            else:
                print("Many amino acids to model on peptide {} ...".format(pep_to_model))
                statsFile.write("The model failed because many aminoacids are required to be predicted\n")
                
    # Close the report file
    statsFile.close()