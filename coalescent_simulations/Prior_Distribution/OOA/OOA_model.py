import msprime
import tskit
import allel; print('scikit-allel', allel.__version__)
import numpy as np
import random
import numpy as np
import sys
import pandas as pd
import os
import matplotlib.pyplot as plt
import stdpopsim

#Generate seeds
import random
seeds = list(range(1,30))

#Set demographic simulation with stdpopsim
species = stdpopsim.get_species("HomSap")
model   = species.get_demographic_model("OutOfAfrica_3G09")
#Create the contig (ie: change the length of the simulation and the recomb rate params)
contig  = species.get_contig("chr22", length_multiplier=0.01)
#Choose sample sizes for each population (YRI/CEU/CHB)
YRI_sample_size = 1000
CEU_sample_size = 1000
CHB_sample_size = 1000
samples = model.get_samples(YRI_sample_size, CEU_sample_size, CHB_sample_size)
engine = stdpopsim.get_engine("msprime")
#Run simulations across replicates
for i in range(1,len(seeds)):
  ts = engine.simulate(model, contig, samples, seed= seeds[i])
  #Measure Fst with tskita
  #ts.Fst([ts.samples(0), ts.samples(1)], indexes = [(0,1)])
  #Convert to haplotype array using scikit allele again noting population structure
  haplotype_array = np.asarray(allel.HaplotypeArray(ts.genotype_matrix()))
  #Split into the different panels
  haplotype_array = np.transpose(haplotype_array)
  #Read out as csv with pandas
  YRI  = pd.DataFrame((haplotype_array[0:YRI_sample_size:,1:300]))
  CEU  = pd.DataFrame((haplotype_array[(YRI_sample_size):(YRI_sample_size+CEU_sample_size):,1:300]))
  CHB  = pd.DataFrame((haplotype_array[(YRI_sample_size+CEU_sample_size):(YRI_sample_size+CEU_sample_size+CHB_sample_size):,1:300]))
  
  #YRI  = pd.DataFrame((haplotype_array))
  #CEU  = pd.DataFrame((haplotype_array))
  #CHB  = pd.DataFrame((haplotype_array))
  
  YRI.to_csv("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/Prior_Distribution/OOA/msprime_data/YRI_replicate_%d.csv" % (i), index=False)    
  CEU.to_csv("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/Prior_Distribution/OOA/msprime_data/CEU_replicate_%d.csv" % (i), index=False)  
  CHB.to_csv("/Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/Prior_Distribution/OOA/msprime_data/CHB_replicate_%d.csv" % (i), index=False)

os.system("scp /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations/Prior_Distribution/OOA/msprime_data/* tutert@rescomp2.well.ox.ac.uk:/well/mcvean/mtutert/thesis/coalescent_coverage/prior_draws/msprime_data/OOA")



