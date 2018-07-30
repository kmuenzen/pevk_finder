import sys
from Bio.Phylo.PAML import codeml

#folder_path = sys.argv[1]
alignment_file = sys.argv[1] # full path
tree_file = sys.argv[2] # full path
m0_out = sys.argv[3] # full output path
estimated_tree_name = sys.argv[4]
final_out = sys.argv[5]

# Run M0 model to get tree
cmlM0 = codeml.Codeml(alignment=alignment_file, tree=tree_file, out_file=m0_out)
cmlM0.set_options(seqtype=1)
cmlM0.set_options(model=0)
cmlM0.set_options(NSsites=[0])
cmlM0.set_options(omega=0.5)
cmlM0.set_options(CodonFreq=2)
cmlM0.set_options(ndata=1)
cmlM0.set_options(fix_alpha=1)
cmlM0.set_options(Small_Diff=5e-7)

# Run the M0 model
cmlM0.run(command="/Users/kmoney/Documents/paml4.9e/bin/codeml")

# Get tree from m0 results
m0result = codeml.read(m0_out)
NSsites_dict = m0result.get("NSsites")
NSsites0_dict = NSsites_dict.get(0)
estimated_tree = NSsites0_dict.get("tree")

# Write tree to output tree file
f = open(estimated_tree_name, "w")
f.write(estimated_tree)
f.close()

# Now run all sites models
cml = codeml.Codeml(alignment=alignment_file, tree=estimated_tree_name, out_file=final_out)
cml.set_options(seqtype=1)
cml.set_options(model=0)
cml.set_options(NSsites=[0,1,2,7,8])
cml.set_options(omega=[0.5,1,3])
cml.set_options(CodonFreq=2)
cml.set_options(ndata=1)
cml.set_options(fix_alpha=1)
cml.set_options(Small_Diff=5e-7)

# Run all sites models
cml.run(command="/Users/kmoney/Documents/paml4.9e/bin/codeml")
