#!/usr/bin/env python

# (c) Copyright Rosetta Commons Member Institutions. 
# (c) This file is part of the Rosetta software suite and is made available under license. 
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons. 
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be 
# (c) addressed to University of Washington CoMotion, email: license@uw.edu. 
"""Brief:   This PyRosetta script does blah.blah blah blah 

Params:  ./blah.py .pdb 

Example: ./blah.py foo.pdb 1000

Remarks: Blah blah blah, blah, blah.

Author:  Jason W. Labonte

"""
## Import Standard Python libraries
import argparse
## Next set of imports, more space efficient 
import pyrosetta 
from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover
from pyrosetta.rosetta.core.pose.carbohydrates import glycosylate_pose
from pyrosetta.rosetta.core.pack.task.operation import InitializeFromCommandline, RestrictToRepacking
from pyrosetta.rosetta.core.select.residue_selector import ResidueIndexSelector, NeighborhoodResidueSelector
from pyrosetta.rosetta.core.simple_metrics.metrics import RMSDMetric
from pyrosetta.rosetta.core.pack.task.operation import OperateOnResidueSubset, PreventRepackingRLT
from pyrosetta.rosetta.protocols.rigid import RigidBodyPerturbNoCenterMover as RBPNCMover 
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.rosetta.core.select.movemap import MoveMapFactory
## Need to do something with them 
from pyrosetta.rosetta.core.scoring.constraints import *
from pyrosetta.rosetta.protocols.carbohydrates import *
from pyrosetta.rosetta.core.scoring import *
#from pyrosetta.rosetta.scoring.constraints import atom_pair_constraint

#from pyrosetta import *
'''
from pyrosetta.rosetta import *
from pyrosetta.rosetta.protocols.carbohydrates import *
from pyrosetta.rosetta.core.pose import *
from pyrosetta.rosetta.core.scoring import *
from pyrosetta.rosetta.protocols.minimization_packing import *
#from pyrosetta.rosetta.core.pose.carbohydrates import glycosylate_pose, glycosylate_pose_by_file
#from pyrosetta.rosetta.protocols.carbohydrates import SimpleGlycosylateMover
from pyrosetta.rosetta.protocols.rosetta_scripts import *
from pyrosetta.rosetta.core.scoring.constraints import *
#from pyrosetta.rosetta.core.select.movemap import MoveMapFactory
'''
## A global dictionary of sugars that I can glycosylate the acceptor peptide with 
sugar_dictionary = {}
sugar_dictionary["core1"] = "b-D-Galp-(1->3)-a-D-GalpNAc"

class OGlycosylationStartingStructureFormation:
    def __init__(self):
        self.acceptor_peptide_sequence = None
        self.acceptor_peptide_sugar_name = None
        self.acceptor_peptide_glycosylation_location = None
        self.reference_pose_file = None 
        self.constraints_file = None
        self.repack_interface = False 
        self.decoy_numbers = None

        # These variables are not user inputs 
        self.acceptor_peptide_pose = None 
        self.input_enzyme_pose = None
        self.output_enzyme_pose = None 
        self.reference_pose = None
        self.sfxn = None
        self.initial_energy = None
        self.acceptor_chain = None
        self.donor_chain = None
        self.sugar_chain_start = None
        self.sugar_chain_end = None
        self.enzyme_chain = None

    def _relax_sugars(self, input_pose):
        sampler = GlycanSampler()
        sampler.set_rounds(100)     
        sampler.apply(input_pose)

    def _prepare_acceptor_peptide_with_glycan(self):
        """ This program takes the sequence of the peptide (the peptide must have one 
        THR residue for glycosylation to happen and should be 2-3 residue long)
        and outputs a pose object of the sequence with the core1 sugar glycosylated 
        on the specified THR residue. """
        self.acceptor_peptide_pose = pyrosetta.pose_from_sequence(self.acceptor_peptide_sequence)
        if self.acceptor_peptide_sugar_name != None: 
            glycosylator = SimpleGlycosylateMover()
        if self.acceptor_peptide_sugar_name not in sugar_dictionary: print("The sugar not known!!")
        glycosylator.set_glycosylation(sugar_dictionary[self.acceptor_peptide_sugar_name])
        glycosylator.set_position(self.acceptor_peptide_glycosylation_location)
        glycosylator.apply(self.acceptor_peptide_pose)
        self._relax_sugars(self.acceptor_peptide_pose) 

    def _merge_two_pose_by_jump(self, enzyme_pose, acceptor_peptide_pose):
        """merge the base pose into enzyme pose without making any covalent bond"""
        anchor_point = int(enzyme_pose.total_residue()) - 2 # just to make sure that the UDP-sugar remains the last residues in the whole pose
        enzyme_pose.append_pose_by_jump(acceptor_peptide_pose, anchor_point)

    def _make_fold_tree_docking_friendly(self, enzyme_pose):
        """
        The function changes the fold_tree of the pose from the default one to a more dockign friendly. 
        """
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
                self.enzyme_chain = i+1
        
        # Criteria 2 and 3: the acceptor peptide will not have any sugar in it and the donor molecule will have UDP it it 
        for i in range(0,len(chain_length_list)):
            if (i+1 != self.enzyme_chain):
                assumption = True 
                for j in range(chain_begin_list[i],chain_end_list[i]+1):
                    name = str(enzyme_pose.residue(j).name())
                    name = name[0:3]
                    if name == 'UDP': self.donor_chain = i+1
                    if name not in AA_set:
                        assumption = False
                if assumption == True: 
                    self.acceptor_chain = i+1
        
        # Finding the pivot residue 
        reference = enzyme_pose.residue(chain_begin_list[self.donor_chain-1]).xyz("3OPB")
        distance_lowest = float('inf')
        pivot_residue = float('inf')
        for i in range(chain_begin_list[self.enzyme_chain-1], chain_end_list[self.enzyme_chain-1]+1):
            residue1 = enzyme_pose.residue(i).xyz("CA")
            distance_curr = (residue1-reference).norm()
            if distance_curr < distance_lowest:
                distance_lowest = distance_curr
                pivot_residue = i

        # Now actually changing the fold tree in an automated fashion 
        first_residue_enzyme = int(enzyme_pose.chain_begin(self.enzyme_chain))
        last_residue_enzyme = int(enzyme_pose.chain_end(self.enzyme_chain))
        
        anchor_residue_enzyme = pivot_residue
        
        first_residue_substrate = int(enzyme_pose.chain_begin(self.acceptor_chain))
        last_residue_substrate = int(enzyme_pose.chain_end(self.acceptor_chain))
        anchor_residue_substrate = int((first_residue_substrate+last_residue_substrate)/2)
        first_donor = int(enzyme_pose.chain_begin(self.donor_chain))
        last_donor = int(enzyme_pose.chain_end(self.donor_chain))

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
            self.sugar_chain_start = int(enzyme_pose.chain_begin(sugar_chain))
            self.sugar_chain_end = int(enzyme_pose.chain_end(sugar_chain))
            ft_docking.add_edge(self.sugar_chain_start,self.sugar_chain_end,-1)

        print(enzyme_pose.fold_tree())
        enzyme_pose.fold_tree(ft_docking)
        print(enzyme_pose.fold_tree())

    def _large_rigid_body_moves(self, curr_enzyme_pose):
        rigid_body = RBPNCMover()
        rigid_body.add_jump(1)
        rigid_body.rot_magnitude(10)
        rigid_body.trans_magnitude(5.0)

        first_pose = curr_enzyme_pose.clone()
        num_iter = 0
        while (self.sfxn(curr_enzyme_pose) > 0 ):
            num_iter += 1
            rigid_body.apply(first_pose)
            if (self.sfxn(first_pose) < self.sfxn(curr_enzyme_pose)):
                curr_enzyme_pose = first_pose.clone()
            else:
                first_pose = curr_enzyme_pose.clone()
            if (num_iter > 500):
                print("Iternation count exceeded !!")
                break
        self.sfxn.show(curr_enzyme_pose)
        return curr_enzyme_pose

    def _finer_rigid_body_move(self, curr_enzyme_pose):
        rigid_body_trans = RBPNCMover()
        rigid_body_trans.add_jump(1)
        rigid_body_trans.rot_magnitude(0)
        rigid_body_trans.trans_magnitude(0.1)
        
        rigid_body_rot = RBPNCMover()
        rigid_body_rot.add_jump(1)
        rigid_body_rot.rot_magnitude(1)
        rigid_body_rot.trans_magnitude(0)
        
        first_pose = curr_enzyme_pose.clone()
        num_iter = 0
        while (self.sfxn(curr_enzyme_pose) > (self.initial_energy) ):
            num_iter += 1
            rigid_body_rot.apply(first_pose)
            if (self.sfxn(first_pose) < self.sfxn(curr_enzyme_pose)):
                curr_enzyme_pose = first_pose.clone()
            rigid_body_trans.apply(first_pose)
            if (self.sfxn(first_pose) < self.sfxn(curr_enzyme_pose)):
                curr_enzyme_pose = first_pose.clone()
            else:
                first_pose = curr_enzyme_pose.clone()
            if (num_iter > 500):
                print("Iternation count exceeded !!")
                break
        self.sfxn.show(curr_enzyme_pose)
        return curr_enzyme_pose

    def _relax_only_acceptor_peptide(self, curr_enzyme_pose):
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(InitializeFromCommandline())
        tf.push_back(RestrictToRepacking())

        acceptor_peptide_selector = ResidueIndexSelector()
        for j in range(int(curr_enzyme_pose.chain_begin(self.acceptor_chain)),int(curr_enzyme_pose.total_residue()) + 1): # Can be generalised 
            acceptor_peptide_selector.append_index(j)
        acceptor_peptide_selector.apply(curr_enzyme_pose)

        prevent_repackign_rlt = PreventRepackingRLT()
        prevent_subset_repacking = OperateOnResidueSubset(prevent_repackign_rlt, acceptor_peptide_selector, flip_subset=True)

        tf.push_back(prevent_subset_repacking)
        mmf = MoveMapFactory()
        mmf.all_bb(False)
        mmf.all_chi(True)
        mmf.add_bb_action(pyrosetta.rosetta.core.select.movemap.move_map_action.mm_enable, acceptor_peptide_selector)
        print(tf.create_task_and_apply_taskoperations(curr_enzyme_pose))

        fr = FastRelax()
        fr.set_task_factory(tf)
        fr.set_scorefxn(self.sfxn)
        fr.set_movemap_factory(mmf)
        fr.max_iter(0)
        fr.apply(curr_enzyme_pose)    
        self.sfxn.show(curr_enzyme_pose)   

    def _relax_the_interface(self, curr_enzyme_pose):
        """This function packs the rotamers in the interface region of the acceptor
        peptide and the enzyme."""
        tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
        tf.push_back(InitializeFromCommandline())
        tf.push_back(RestrictToRepacking())

        acceptor_peptide_selector = ResidueIndexSelector()
        acceptor_peptide_neighborhood_selector = NeighborhoodResidueSelector()
        for j in range(int(curr_enzyme_pose.chain_begin(self.acceptor_chain)),int(curr_enzyme_pose.total_residue()) + 1): # Can be generalised 
            acceptor_peptide_selector.append_index(j)
        acceptor_peptide_selector.apply(curr_enzyme_pose)
        acceptor_peptide_neighborhood_selector.set_focus_selector(acceptor_peptide_selector)
        acceptor_peptide_neighborhood_selector.set_include_focus_in_subset(True)

        prevent_repackign_rlt = PreventRepackingRLT()
        prevent_subset_repacking = OperateOnResidueSubset(prevent_repackign_rlt, acceptor_peptide_neighborhood_selector, flip_subset=True)

        tf.push_back(prevent_subset_repacking)
        mmf = MoveMapFactory()
        mmf.all_bb(False)
        mmf.all_chi(True)
        mmf.add_bb_action(pyrosetta.rosetta.core.select.movemap.move_map_action.mm_enable, acceptor_peptide_neighborhood_selector)
        print(tf.create_task_and_apply_taskoperations(curr_enzyme_pose))

        fr = FastRelax()
        fr.set_task_factory(tf)
        fr.set_scorefxn(self.sfxn)
        fr.set_movemap_factory(mmf)
        fr.max_iter(0)
        fr.apply(curr_enzyme_pose)    
        self.sfxn.show(curr_enzyme_pose)   

    def _add_acceptor_peptide_and_enzyme(self, decoy_numbers):
        """ This program takes the base sugar pose and the pose of the enzyme with donor moiety
        and the constraints and make sure the they obey the constraint in subsequent relax
        procedure as a procedure to bypass manual overlaying using PyMOL. """
        
        self.sfxn = pyrosetta.create_score_function("ref2015_cart.wts")
        self.initial_energy = self.sfxn(self.input_enzyme_pose)

        # merging two poses to form one 
        self._merge_two_pose_by_jump(self.input_enzyme_pose, self.acceptor_peptide_pose)
        
        # Adding the constraints to the score term
        self.sfxn.set_weight(atom_pair_constraint, 1.0)
        add_fa_constraints_from_cmdline(self.input_enzyme_pose, self.sfxn)
        print(self.sfxn.show(self.input_enzyme_pose))

        # Change the FoldTree to suite Docking protocol 
        self._make_fold_tree_docking_friendly(self.input_enzyme_pose)

        # Defining the reference_pose for peptide-RMSD calculation
        # and Initializing the output lists
        if self.reference_pose_file:
            self.reference_pose = pyrosetta.pose_from_file(self.reference_pose_file)
            peptide_reference_pose = self.input_enzyme_pose.clone()
            score_list = []  
            distance_list = []
            core1_rmsd_list = []
            peptide_rmsd_list = []

        minimum_score = float('inf')
        self.output_enzyme_pose = self.input_enzyme_pose.clone()
        for trial_number in range(0,decoy_numbers):
            print("Decoy Number: " + str(trial_number))
            curr_enzyme_pose = self.input_enzyme_pose.clone()
            curr_enzyme_pose = self._large_rigid_body_moves(curr_enzyme_pose)
            curr_enzyme_pose = self._finer_rigid_body_move(curr_enzyme_pose)
            if self.repack_interface not in [True, False]:
                print("repack_interface must be a boolean True or False !!!")
            if self.repack_interface == False:
                self._relax_only_acceptor_peptide(curr_enzyme_pose)
            else:
                self._relax_the_interface(curr_enzyme_pose)
            # take the acceptor_peptide_with_sugar 10 A back --------------------------------------------------------------
            '''
            translation_magnitude = 10
            rigid_body_jump = 1

            translate_away = pyrosetta.rosetta.protocols.rigid.RigidBodyTransMover(curr_enzyme_pose, rigid_body_jump)
            translate_away.step_size(translation_magnitude)
            translate_away.apply(curr_enzyme_pose)
            # dock_perturb the acceptor_peptide_with_sugar ----------------------------------------------------------------
            
            # slide it into contact ---------------------------------------------------------------------------------------
            slide = pyrosetta.rosetta.protocols.dokcing.FaDockingSlideIntoContact(rigid_body_jump, slide_axis)
            slide.apply(curr_enzyme_pose)
            '''
            # --------------------------------------------------------------------------------------------------------------
            if self.sfxn(curr_enzyme_pose) < minimum_score:
                self.output_enzyme_pose = curr_enzyme_pose.clone()
                minimum_score = self.sfxn(curr_enzyme_pose)

            ## RMSD and stuff like that for benchmarking purpose, a standard user need not bother about this part 
            if self.reference_pose_file:
                score_list.append(self.sfxn(curr_enzyme_pose))
                atom1 = curr_enzyme_pose.residue(int(curr_enzyme_pose.chain_end(self.donor_chain))).xyz("C1") 
                atom2 = curr_enzyme_pose.residue(int(curr_enzyme_pose.chain_end(self.acceptor_chain))+1).xyz("O6") 
                distance_list.append((atom1-atom2).norm())
                ## RMSD of the core 1 sugar -------------------------------------------------------------------- 
                allEnzyme_experimental = ResidueIndexSelector()
                allEnzyme_reference = ResidueIndexSelector()
                for j in range(int(curr_enzyme_pose.chain_begin(self.enzyme_chain)), int(curr_enzyme_pose.chain_end(self.enzyme_chain))+1):
                    allEnzyme_experimental.append_index(j)
                    allEnzyme_reference.append_index(j)

                core1_experimental = ResidueIndexSelector()
                for j in range(self.sugar_chain_start,self.sugar_chain_end+1):
                    core1_experimental.append_index(j)
                core1_reference = ResidueIndexSelector()
                for j in range(372,374): # Do i need to generalize it??? 
                    core1_reference.append_index(j)

                core1_rmsd = RMSDMetric(self.reference_pose)
                core1_rmsd.set_residue_selector(core1_experimental)
                core1_rmsd.set_residue_selector_reference(core1_reference)
                core1_rmsd.set_residue_selector_super(allEnzyme_experimental)
                core1_rmsd.set_residue_selector_super_reference(allEnzyme_reference)
                core1_rmsd.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_all_heavy)
                core1_rmsd_list.append(float(core1_rmsd.calculate(curr_enzyme_pose)))
                ## ------------------------------------------------------------------------------------------------
                ## RMSD of the peptide final pose with its initial form
                peptide_experimental = ResidueIndexSelector()
                peptide_reference = ResidueIndexSelector()
                for j in range(int(curr_enzyme_pose.chain_begin(self.acceptor_chain)), int(curr_enzyme_pose.chain_end(self.acceptor_chain))+1):
                    peptide_experimental.append_index(j)
                    peptide_reference.append_index(j)
                
                peptide_rmsd = RMSDMetric(peptide_reference_pose)
                peptide_rmsd.set_residue_selector(peptide_experimental)
                peptide_rmsd.set_residue_selector_reference(peptide_reference)
                peptide_rmsd.set_residue_selector_super(allEnzyme_experimental)
                peptide_rmsd.set_residue_selector_super_reference(allEnzyme_reference)
                peptide_rmsd.set_rmsd_type(pyrosetta.rosetta.core.scoring.rmsd_all)
                peptide_rmsd_list.append(float(peptide_rmsd.calculate(curr_enzyme_pose)))
                ## ---------------------------------------------------------------------------------------------------

        ## Printing out the benchmarking inforfmation and storing it inside a file for further analysis, if required 
        print(self.sfxn.show(self.output_enzyme_pose))
        if self.reference_pose_file:
            print("score_list = ", score_list)
            print("distance_list = ", distance_list)
            print("core1_rmsd_list = ", core1_rmsd_list)
            print("peptide_rmsd_list = ", peptide_rmsd_list)

            output_file_name = "O_glycosylation_starting_structure_benchmark_" + self.acceptor_peptide_sequence + ".txt"
            file_opened = open(output_file_name, 'w')
            str1 = "score_list = " + str(score_list) + "\n"
            str2 = "distance_list = " + str(distance_list) + "\n"
            str3 = "core1_rmsd_list = " + str(core1_rmsd_list) + "\n"
            str4 = "peptide_rmsd_list = " + str(peptide_rmsd_list) + "\n"
            file_opened.write(str1)
            file_opened.write(str2)
            file_opened.write(str3)
            file_opened.write(str4)
            file_opened.close()

    def apply(self, input_pose):
        """
        This is the function the user should use to apply the full "Mover" on the input_pose which 
        is actyally the enzyme + donor pose and return the pose with enzyme + donor residue + 
        acceptor residue.
        """
        self.input_enzyme_pose = input_pose
        self._prepare_acceptor_peptide_with_glycan()
        self._add_acceptor_peptide_and_enzyme(self.decoy_numbers)

def argument_parsing():
    info="""
        Will write something here after some time
    """
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-pepseq','--peptide_sequence', type=str, required=True, 
        help="Input the one letter amino acid code of the acceptor peptide sequence.")
    parser.add_argument('-sugar', '--sugar_name', type=str, required=True,
        help="Input the name of the sugar you want to use, e.g. core1.")
    parser.add_argument('-gly_loc', '--glycosylation_location', type=int, required=True,
        help="Input the index (starting from 1) of the amino acid in the acceptor \
        peptide sequence that will be glycosylated by the input sugar")
    parser.add_argument('-ref', '--reference_file', type=str, required=False,
        help="If you want the Benchmarking to be done, input the location of \
        the reference pdb file.")
    parser.add_argument('-rI', '--repack_interface', type=bool, required=False, 
        default=False, help="Input if you want to repack the interface region \
        of the acceptor peptide and the enzyme.")
    parser.add_argument('-enzyme', '--enzyme_file', type=str, required=True, 
        help="Input the location of the enzyme + donor pdb file.")
    parser.add_argument('-cnstr', '--constraints_file', type=str, required=True,
        help="Input the lcoation of the constraints file.")
    parser.add_argument('-decoy', '--decoy_numbers', type=int, required=False,
        default=1, help="Input the number of independent decoys the user wants.")

    return parser.parse_args()

if __name__=="__main__":
    args = argument_parsing()

    mover = OGlycosylationStartingStructureFormation()
    mover.acceptor_peptide_sequence = args.peptide_sequence
    mover.acceptor_peptide_sugar_name = args.sugar_name
    mover.acceptor_peptide_glycosylation_location = args.glycosylation_location 
    mover.reference_pose_file = args.reference_file
    mover.constraints_file = args.constraints_file
    mover.repack_interface = args.repack_interface
    mover.decoy_numbers = args.decoy_numbers

    init_flags =  "-include_sugars -maintain_links -auto_detect_glycan_connections -alternate_3_letter_codes pdb_sugar" + " -constraints:cst_fa_file " + mover.constraints_file
    pyrosetta.init(init_flags)

    enzyme_pose = pyrosetta.pose_from_pdb(args.enzyme_file)
    mover.apply(enzyme_pose)


"""
if __name__=="__main__":
    mover = OGlycosylationStartingStructureFormation()
    mover.acceptor_peptide_sequence = "ETTSHST"
    mover.acceptor_peptide_sugar_name = "core1"
    mover.acceptor_peptide_glycosylation_location = 3
    mover.reference_pose_file = "/home/souvadra/myGitFolders/Glycosyltransferase/Acceptor-Donor-Enzyme/GlcNAc-added-before-GalBGalNAc/3OTK-closed-monomer-alpha-GlcNAc_2GAM-GalBGalNAc.pdb"
    mover.constraints_file = "3OTK_constraints_file.cst"
    mover.repack_interface = True 
    mover.decoy_numbers = 1

    init_flags =  "-include_sugars -maintain_links -auto_detect_glycan_connections -alternate_3_letter_codes pdb_sugar" + " -constraints:cst_fa_file " + mover.constraints_file
    pyrosetta.init(init_flags)
    enzyme_pose_file = "/home/souvadra/myGitFolders/Glycosyltransferase/glycan_sampler_pipeline/output_3OTK-closed-S217C/3OTK-closed-monomer-alpha-GlcNAc-S217C_0005_598_manually_removed.pdb"
    enzyme_pose = pyrosetta.pose_from_pdb(enzyme_pose_file)
    mover.apply(enzyme_pose)

    #mover._prepare_acceptor_peptide_with_glycan()
    #mover._add_acceptor_peptide_and_enzyme(2)
"""