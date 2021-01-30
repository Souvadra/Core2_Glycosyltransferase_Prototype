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

from pyrosetta.rosetta.core.scoring.constraints import *
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory

def mergePoses(base_pose, enzyme_pose):
    """merge the base pose into enzyme pose without making any covalent bond"""
    anchor_point = int(enzyme_pose.total_residue()) - 2 # just to make sure that the UDP-sugar remains the last residues in the whole pose
    enzyme_pose.append_pose_by_jump(base_pose, anchor_point)

def addBaseSugarAndEnzyme(base_pose, enzyme_pose, constraints_file, decoy_numbers=5):
    """ This program takes the base sugar pose and the pose of the enzyme with donor moiety
        and the constraints and make sure the they obey the constraint in subsequent relax
        procedure as a procedure to bypass manual overlaying using PyMOL. """
    init("-constraints:cst_fa_file " + constraints_file)
    
    sfxn = pyrosetta.create_score_function("ref2015_cart.wts") #get_fa_scorefxn()
    initial_energy = sfxn(enzyme_pose)
    # merging two poses to form one 
    mergePoses(base_pose, enzyme_pose)

    sfxn.set_weight(atom_pair_constraint, 1.0)
    add_fa_constraints_from_cmdline(enzyme_pose, sfxn)
    print(sfxn.show(enzyme_pose))

    ## Change the FoldTree to more Docking friendly option as suggested by Pooja 
    first_residue_enzyme = 1
    last_residue_enzyme = 369
    anchor_residue_enzyme = 130
    anchor_residue_substrate = 373
    first_residue_substrate = 372
    last_residue_substrate = 374
    first_donor = 370
    last_donor = 371

    ft_docking = pyrosetta.FoldTree()
    ft_docking.add_edge(anchor_residue_enzyme,1,-1)
    ft_docking.add_edge(anchor_residue_enzyme,last_residue_enzyme,-1)

    ft_docking.add_edge(anchor_residue_enzyme, anchor_residue_substrate, 1) 
    ft_docking.add_edge(anchor_residue_enzyme, first_donor, 2)
            
    ft_docking.add_edge(anchor_residue_substrate, first_residue_substrate, -1) # downstream base peptide 
    ft_docking.add_edge(anchor_residue_substrate, last_residue_substrate, -1) # upstream base peptide 
    ft_docking.add_edge(first_donor,last_donor,-1) # the donor UDP-GlcNAc backbone

    chemical_edges = enzyme_pose.fold_tree().get_chemical_edges() # connections to the glycan residues from anchor out 
    for chemical in range(1, len(chemical_edges) + 1):
        ft_docking.add_edge(chemical_edges[chemical].start(), chemical_edges[chemical].stop(), -2)

    ft_docking.add_edge(375,376,-1) # Adding the covalent bond between the two sugars of core 1 
    print(enzyme_pose.fold_tree())
    enzyme_pose.fold_tree(ft_docking)
    print(enzyme_pose.fold_tree())
    
    ## setting the movemap 
    mm = pyrosetta.rosetta.core.kinematics.MoveMap()
    mm.set_bb(False)
    mm.set_chi(False)
    mm.set_jump(True)

    RMSD_list = []
    score_list = []  
    referencePose = pose_from_pdb("/home/souvadra/myGitFolders/Glycosyltransferase/Acceptor-Donor-Enzyme/GlcNAc-added-before-GalBGalNAc/3OTK-closed-monomer-alpha-GlcNAc_2GAM-GalBGalNAc.pdb")
    for trial_number in range(0,decoy_numbers):
        print("Decoy number: " + str(trial_number))
        ## Use Rigid Body Mover to bring the base peptide close to the enzyme
        curr_enzyme_pose = enzyme_pose.clone()
        rigid_body = pyrosetta.rosetta.protocols.rigid.RigidBodyPerturbNoCenterMover()
        rigid_body.add_jump(1)
        rigid_body.rot_magnitude(10)
        rigid_body.trans_magnitude(5.0)
        
        first_pose = curr_enzyme_pose.clone()
        num_iter = 0
        while (sfxn(curr_enzyme_pose) > 0 ):
            num_iter += 1
            rigid_body.apply(first_pose)
            if (sfxn(first_pose) < sfxn(curr_enzyme_pose)):
                curr_enzyme_pose = first_pose.clone()
                #print(sfxn(curr_enzyme_pose))
            else:
                first_pose = curr_enzyme_pose.clone()
            if (num_iter > 1000):
                print("Iternation count exceeded !!")
                break
            #print(sfxn(enzyme_pose))
        sfxn.show(curr_enzyme_pose)
        #print(curr_enzyme_pose.fold_tree())

        rigid_body2 = pyrosetta.rosetta.protocols.rigid.RigidBodyPerturbNoCenterMover()
        rigid_body2.add_jump(1)
        rigid_body2.rot_magnitude(0)
        rigid_body2.trans_magnitude(0.1)
        
        rigid_body_rot = pyrosetta.rosetta.protocols.rigid.RigidBodyPerturbNoCenterMover()
        rigid_body_rot.add_jump(1)
        rigid_body_rot.rot_magnitude(1)
        rigid_body_rot.trans_magnitude(0)
        
        first_pose = curr_enzyme_pose.clone()
        num_iter = 0
        while (sfxn(curr_enzyme_pose) > (initial_energy) ):
            num_iter += 1
            rigid_body_rot.apply(first_pose)
            if (sfxn(first_pose) < sfxn(curr_enzyme_pose)):
                curr_enzyme_pose = first_pose.clone()
                #print(sfxn(curr_enzyme_pose), "  rotation-iteration ")    
            rigid_body2.apply(first_pose)
            if (sfxn(first_pose) < sfxn(curr_enzyme_pose)):
                curr_enzyme_pose = first_pose.clone()
                #print(sfxn(curr_enzyme_pose), "  translation-iteration ")
            else:
                first_pose = curr_enzyme_pose.clone()
            if (num_iter > 500):
                print("Iternation count exceeded !!")
                break
            #print(sfxn(enzyme_pose))
        score_list.append(sfxn(curr_enzyme_pose))
        RMSD_list.append(CA_rmsd(curr_enzyme_pose,referencePose))
        dumping_name = "merging_result" + str(trial_number) + ".pdb"
        curr_enzyme_pose.dump_pdb(dumping_name)

    '''
    ## Apply the minimizer (use Cartesian coordinates)
    minimizer =  pyrosetta.rosetta.protocols.minimization_packing.MinMover()
    minimizer.cartesian(True)
    minimizer.tolerance(0.001) 
    minimizer.movemap(mm) 
    minimizer.min_type("lbfgs_armijo_nonmonotone")
    minimizer.score_function(sfxn)
    minimizer.apply(enzyme_pose)
    '''
    min_energy = min(score_list)
    for i in range(0,len(score_list)):
        if score_list[i] == min_energy:
            print("merging_result" + str(i) + ".pdb has the lowest energy.")
            output_pose = pose_from_pdb("merging_result" + str(i) + ".pdb")

    print(sfxn.show(output_pose))
    output_pose.dump_pdb("merging_result2.pdb")
    #return None # So far, will change this shortly
    return output_pose