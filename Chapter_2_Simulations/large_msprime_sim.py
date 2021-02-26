#Run a giant msprime simulation with lots and lots of mutations and 'two' populations, both of sample size 1000
import msprime
import allel; print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd

#Run the msprime simulation
map_positions = [i*1 for i in range(0, 3)]
map_rates = [0,1e-1,0]
my_map = msprime.RecombinationMap(map_positions, map_rates)
tree = msprime.simulate(sample_size=2000, Ne=1000,recombination_map = my_map, mutation_rate=5e-1)
haps = np.transpose(np.asarray(allel.HaplotypeArray(tree.genotype_matrix())))
np.shape(haps)
#Write out the data
#Extract the GWAS file
pop_1 = pd.DataFrame(haps[0:1000,])
pop_2 = pd.DataFrame(haps[1001:2000,])
pop_1.to_csv("msprime_data/conditional_ld/pop_1.csv", index = False)
pop_2.to_csv("msprime_data/conditional_ld/pop_2.csv", index = False)
