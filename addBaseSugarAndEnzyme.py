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
    chain_begin_list = []
    chain_end_list = []
    chain_begin_list.append(int(enzyme_pose.chain_begin(1)))
    chain_end_list.append(int(enzyme_pose.chain_end(1)))
    while (chain_end_list[-1] < int(enzyme_pose.total_residue())):
        chain_num = int(enzyme_pose.chain(chain_end_list[-1]+1))
        chain_begin_list.append(int(enzyme_pose.chain_begin(chain_num)))
        chain_end_list.append(int(enzyme_pose.chain_end(chain_num)))
    total_chains = len(chain_begin_list)

    AA_set = {'ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'}
    
    # Criteria 1: Enzyme generally the largest chain of this reaction
    chain_length_list = []
    for i in range(0,len(chain_begin_list)):
        chain_length_list.append(chain_end_list[i]-chain_begin_list[i])
    largest_chain_length = max(chain_length_list)

    for i in range(0,len(chain_length_list)):
        if chain_length_list[i] == largest_chain_length:
            enzyme_chain = i+1
    
    # Criteria 2 and 3: the acceptor peptide will not have any sugar in it and the donor molecule will have UDP it it 
    for i in range(0,len(chain_length_list)):
        if (i+1 != enzyme_chain):
            assumption = True 
            for j in range(chain_begin_list[i],chain_end_list[i]+1):
                name = str(enzyme_pose.residue(j).name())
                name = name[0:3]
                if name == 'UDP': donor_chain = i+1
                if name not in AA_set:
                    assumption = False
            if assumption == True: 
                acceptor_chain = i+1
    
    # Finding the pivot residue 
    reference = enzyme_pose.residue(chain_begin_list[donor_chain-1]).xyz("3OPB")
    distance_lowest = float('inf')
    pivot_residue = float('inf')
    for i in range(chain_begin_list[enzyme_chain-1], chain_end_list[enzyme_chain-1]+1):
        residue1 = enzyme_pose.residue(i).xyz("CA")
        distance_curr = (residue1-reference).norm()
        if distance_curr < distance_lowest:
            distance_lowest = distance_curr
            pivot_residue = i

    # Now actually changing the fold tree in an automated fashion 
    first_residue_enzyme = int(enzyme_pose.chain_begin(enzyme_chain))
    last_residue_enzyme = int(enzyme_pose.chain_end(enzyme_chain))
    
    anchor_residue_enzyme = pivot_residue
    
    first_residue_substrate = int(enzyme_pose.chain_begin(acceptor_chain))
    last_residue_substrate = int(enzyme_pose.chain_end(acceptor_chain))
    anchor_residue_substrate = int((first_residue_substrate+last_residue_substrate)/2)
    first_donor = int(enzyme_pose.chain_begin(donor_chain))
    last_donor = int(enzyme_pose.chain_end(donor_chain))

    ft_docking = pyrosetta.FoldTree()
    ft_docking.add_edge(anchor_residue_enzyme,first_residue_enzyme,-1)
    ft_docking.add_edge(anchor_residue_enzyme,last_residue_enzyme,-1)

    ft_docking.add_edge(anchor_residue_enzyme, anchor_residue_substrate, 1) 
    ft_docking.add_edge(anchor_residue_enzyme, first_donor, 2)
            
    ft_docking.add_edge(anchor_residue_substrate, first_residue_substrate, -1) # downstream base peptide 
    ft_docking.add_edge(anchor_residue_substrate, last_residue_substrate, -1) # upstream base peptide 
    ft_docking.add_edge(first_donor,last_donor,-1) # the donor UDP-GlcNAc backbone

    chemical_edges = enzyme_pose.fold_tree().get_chemical_edges() # connections to the glycan residues from anchor out 
    for chemical in range(1, len(chemical_edges) + 1):
        ft_docking.add_edge(chemical_edges[chemical].start(), chemical_edges[chemical].stop(), -2)
        sugar_start_res = chemical_edges[chemical].stop()
        sugar_chain = int(enzyme_pose.chain(sugar_start_res))
        sugar_chain_start = int(enzyme_pose.chain_begin(sugar_chain))
        sugar_chain_end = int(enzyme_pose.chain_end(sugar_chain))
        ft_docking.add_edge(sugar_chain_start,sugar_chain_end,-1)

    # ft_docking.add_edge(375,376,-1) # Adding the covalent bond between the two sugars of core 1 
    print(enzyme_pose.fold_tree())
    enzyme_pose.fold_tree(ft_docking)
    print(enzyme_pose.fold_tree())
    
    score_list = []  
    distance_list = []
    minimum_score = float('inf')
    answer_pose = enzyme_pose.clone()
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
            if (num_iter > 500):
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
        
        ## Final Relax protocol, to make sure things get proper
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
        tf.push_back(pyrosetta.rosetta.core.pack.task.operation.RestrictToRepacking())

        #packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
        #packer.task_factory(tf)

        acceptor_peptide_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        for j in range(int(curr_enzyme_pose.chain_begin(acceptor_chain)),int(enzyme_pose.total_residue()) + 1): # Can be generalised 
            acceptor_peptide_selector.append_index(j)
        acceptor_peptide_selector.apply(curr_enzyme_pose)

        prevent_repackign_rlt = pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT()
        prevent_subset_repacking = pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(prevent_repackign_rlt, acceptor_peptide_selector, flip_subset=True)

        tf.push_back(prevent_subset_repacking)
        mmf = MoveMapFactory()
        mmf.all_bb(False)
        mmf.all_chi(True)
        mmf.add_bb_action(pyrosetta.rosetta.core.select.movemap.move_map_action.mm_enable, acceptor_peptide_selector)

        fr = pyrosetta.rosetta.protocols.relax.FastRelax()
        fr.set_task_factory(tf)
        fr.set_scorefxn(sfxn)
        fr.set_movemap_factory(mmf)
        fr.max_iter(0)
        fr.apply(curr_enzyme_pose)
        ## Relax protocol ends 
        
        if sfxn(curr_enzyme_pose) < minimum_score:
            answer_pose = curr_enzyme_pose.clone()
            minimum_score = sfxn(curr_enzyme_pose)

        score_list.append(sfxn(curr_enzyme_pose))
        atom1 = curr_enzyme_pose.residue(int(enzyme_pose.chain_end(donor_chain))).xyz("C1") # Can be generalised 
        atom2 = curr_enzyme_pose.residue(int(enzyme_pose.chain_end(acceptor_chain))+1).xyz("O6") # Can be generalised 
        distance_list.append((atom1-atom2).norm())
        #dumping_name = "merging_result" + str(trial_number) + ".pdb"
        #curr_enzyme_pose.dump_pdb(dumping_name)

    print(sfxn.show(answer_pose))
    print("score_list: ", score_list)
    print("distance_list: ", distance_list)
    
    return answer_pose
