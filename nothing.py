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

pyrosetta.init('-include_sugars -write_pdb_link_records -ex1 -ex2 -alternate_3_letter_codes pdb_sugar -auto_detect_glycan_connections -min_bond_length 1.1 -max_bond_length 1.5 -ignore_zero_occupancy false -ignore_unrecognized_res false')
sfxn = pyrosetta.create_score_function("ref2015_cart.wts")
pose = pose_from_pdb("merging_result1.pdb")
sfxn.show(pose)
