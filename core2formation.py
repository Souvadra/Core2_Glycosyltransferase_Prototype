## Rosetta Imports
import sys 
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.carbohydrates import *
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.minimization_packing import *
from pyrosetta.rosetta.core.pose.carbohydrates import glycosylate_pose, glycosylate_pose_by_file
from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover
from pyrosetta.rosetta.protocols.rosetta_scripts import *

## Personal Imports
import basePeptidePrep
import addBaseSugarAndEnzyme

## Pyrosetta Initialization
pyrosetta.init('-include_sugars -write_pdb_link_records -ex1 -ex2 -alternate_3_letter_codes pdb_sugar -auto_detect_glycan_connections -min_bond_length 1.1 -max_bond_length 1.5 -ignore_zero_occupancy false -ignore_unrecognized_res false')

## Command Line Arguments 
'''
Command Line Arguments:
# For the base peptide 
- Sequence of the base peptide (2-3 AA long peptides are ideal) : base_seq
- Position of the sequence containing THR that will be glycosylated : base_pos
- the sugar that will be glycosylated on the base peptide : base_sugar

# input pdb files
- The enzyme with donor peptide 

- The constraint file

'''
base_seq = str(sys.argv[1]) #"STP"
base_position = int(sys.argv[2]) #2
base_sugar = "core1"
input_enzyme_file = "/home/souvadra/myGitFolders/Glycosyltransferase/glycan_sampler_pipeline/output_3OTK-closed-S217C/3OTK-closed-monomer-alpha-GlcNAc-S217C_0005.pdb"

reference_pose_file  = "/home/souvadra/myGitFolders/Glycosyltransferase/Acceptor-Donor-Enzyme/GlcNAc-added-before-GalBGalNAc/3OTK-closed-monomer-alpha-GlcNAc_2GAM-GalBGalNAc.pdb"
constraints_file = "3OTK_constraints_file.cst" # "constraints_file.cst"
#constraints_file = "trialConstraint.cst"

## Prepare the base peptide 
base_pose = basePeptidePrep.basePeptideBuild(base_seq, base_position, base_sugar)
print(base_pose)

## Addition of peptide and sugar motif 
enzyme_pose = pose_from_pdb(input_enzyme_file)
sugars_and_enzyme_pose = addBaseSugarAndEnzyme.addBaseSugarAndEnzyme(base_pose, enzyme_pose, constraints_file,25, reference_pose_file, base_seq)
output_name = "deployment_" + "3OTK_" + base_seq + ".pdb"
sugars_and_enzyme_pose.dump_pdb(output_name)

##########################################################################################
#############                    Calculating RMSD                           ##############
##########################################################################################
allEnzyme_experimental = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
allEnzyme_reference = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
for j in range(1,372):
    allEnzyme_experimental.append_index(j)
    allEnzyme_reference.append_index(j)

core1_experimental = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
for j in range(375,377):
    core1_experimental.append_index(j)
core1_reference = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
for j in range(372,374):
    core1_reference.append_index(j)

## Defining the reference pose 
reference_pose = pose_from_pdb(reference_pose_file)

core1_rmsd = pyrosetta.rosetta.core.simple_metrics.metrics.RMSDMetric(reference_pose)
core1_rmsd.set_residue_selector(core1_experimental)
core1_rmsd.set_residue_selector_reference(core1_reference)
core1_rmsd.set_residue_selector_super(allEnzyme_experimental)
core1_rmsd.set_residue_selector_super_reference(allEnzyme_reference)
core1_rmsd.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_all)
my_RMSD = core1_rmsd.calculate(sugars_and_enzyme_pose)
print("The RMSD of the core1 sugar: ", my_RMSD)
