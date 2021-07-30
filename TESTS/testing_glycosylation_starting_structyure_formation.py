import sys 
sys.path.insert(0,'/home/shati/Core2_Glycosyltransferase_Prototype/')

import O_glycosylation_starting_structure_formation as starting
import pyrosetta 
#from pyrosetta.rosetta.core.scoring.constraints import *
#from pyrosetta.rosetta.core.scoring import *

if __name__=="__main__":
    peptide_sequence_list = ['ETTSHST']
    glycosylation_location_list = [3]
    constraints_file_list = ["3OTK_constraints_file_ETTSHST.cst"]
    repack_interface_list = [False]
    enzyme_pose_file_list = ["/home/shati/Glycosyltransferase/glycan_sampler_pipeline/output_3OTK-closed-S217C/3OTK-closed-monomer-alpha-GlcNAc-S217C_0005_598_manually_removed.pdb"]

    mover = starting.OGlycosylationStartingStructureFormation()
    mover.acceptor_peptide_sequence = peptide_sequence_list[0]
    mover.acceptor_peptide_sugar_name =  "core1"
    mover.acceptor_peptide_glycosylation_location = glycosylation_location_list[0]
    mover.reference_pose_file = "/home/shati/Glycosyltransferase/Acceptor-Donor-Enzyme/GlcNAc-added-before-GalBGalNAc/3OTK-closed-monomer-alpha-GlcNAc_2GAM-GalBGalNAc.pdb"
    mover.constraints_file  = constraints_file_list[0]
    mover.repack_interface = repack_interface_list[0]
    mover.decoy_numbers = 1
    mover.output_pdb = True 
    

    init_flags =  "-include_sugars -maintain_links -auto_detect_glycan_connections -alternate_3_letter_codes pdb_sugar" + " -constraints:cst_fa_file " + mover.constraints_file
    pyrosetta.init(init_flags)

    mover.apply(pyrosetta.pose_from_pdb(enzyme_pose_file_list[0]))

