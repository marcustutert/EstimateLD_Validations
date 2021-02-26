import msprime
import tskit
import allel; print('scikit-allel', allel.__version__)
import random
import numpy as np
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
os.system('ls -l')

sys.path.append("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations")
import population_structure 
tree = msprime.simulate(5).first()
print(tree.draw(format="unicode"))
nreplicates     = 100
divergence_time = 0
x = ancestral_pop_split(divergence_time  = divergence_time,
                        seed             = 1,
                        ref_pop_size     = 1000,
                        gwas_pop_size    = 5000,
                        mutation_rate    = 1e-7,
                        length           = 2e4,
                        Ne_A             = 5000,
                        num_replicates   = nreplicates,
                        statistics       = 1)
# fst_value = []                        
# for i in range(0,nreplicates):
#     fst = float(x[i].Fst([x[i].samples(population=1), x[i].samples(population=2)], indexes = [(0,1)]))
#     fst_value.append(fst)
#Store replicates
for i in range(0,nreplicates):
    gwas = pd.DataFrame(x[i][0])
    ref  = pd.DataFrame(x[i][1])

    gwas.to_csv("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/panel_data/GWAS_panel_replicate_%d_split_%d.csv" % (i,divergence_time), index=False)     
    ref.to_csv("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/panel_data/Ref_panel_replicate_%d_split_%d.csv" % (i,divergence_time), index=False)    
    
os.system("scp /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/panel_data/* tutert@rescomp2.well.ox.ac.uk:/well/mcvean/mtutert/thesis/coalescent_coverage/prior_draws/msprime_data/population_split")

