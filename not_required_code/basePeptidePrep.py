## PyRosetta Imports
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.carbohydrates import *
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.minimization_packing import *
from pyrosetta.rosetta.core.pose.carbohydrates import glycosylate_pose, glycosylate_pose_by_file
from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover
from pyrosetta.rosetta.protocols.rosetta_scripts import *

## storing the sugars that my program is compatible with in the form of a dictionary for further 
## expandability
sugar_dicionary = {}
sugar_dicionary["core1"] = "b-D-Galp-(1->3)-a-D-GalpNAc"

def relax_sugars(pose):
    sampler = GlycanSampler()
    sampler.set_rounds(100)
    sampler.apply(pose)

def basePeptideBuild(sequence, res_number, sugar):
    ''' This program takes the sequence of the peptide (the peptide must have one 
        THR residue for glycosylation to happen and should be 2-3 residue long)
        and outputs a pose object of the sequence with the core1 sugar glycosylated 
        on the specified THR residue. '''
    pose1 = pose_from_sequence(sequence)
    if sugar != None:
        glycosylator = SimpleGlycosylateMover()
        if sugar not in sugar_dicionary: print("The sugar not known!!")
        glycosylator.set_glycosylation(sugar_dicionary[sugar])
        glycosylator.set_position(res_number)
        glycosylator.apply(pose1)
        relax_sugars(pose1) 
    return pose1
