filename = "TESTING_acceptor_peptide_pose.pdb"

with open(filename) as f:
    content = f.readlines()

for line in content:
    print(line)
