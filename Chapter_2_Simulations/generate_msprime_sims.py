import msprime
import allel; print('scikit-allel', allel.__version__)
import numpy as np
import pandas as pd

#Call snakemake params
replicates = snakemake.params
print(replicates)
#We need to generate msprime_sims_across_
map_positions = [i*1 for i in range(0, 3)]
map_rates = [0,1e-1,0]
my_map = msprime.RecombinationMap(map_positions, map_rates)
tree = msprime.simulate(sample_size=1000000, Ne=1000,recombination_map = my_map, mutation_rate=5e-2)
haps = np.transpose(np.asarray(allel.HaplotypeArray(tree.genotype_matrix())))

#Extract the GWAS file
gwas = pd.DataFrame(haps[0:1000,])
gwas.to_csv("msprime_data/gwas_test.csv", index = False)

#Split up the other files into 1000 sample based chunks
for i in range(0,99):
    # TODO: write code...
    ref  = pd.DataFrame(haps[i*1000:((i+1)*1000),])
    ref.to_csv("msprime_data/ref_%s.csv" % i, index = False)
