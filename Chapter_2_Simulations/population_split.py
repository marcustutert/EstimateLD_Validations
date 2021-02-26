#Functions that will generate my ancestral population split
import numpy as np
import allel
import math as math
import msprime

def ancestral_pop_split(divergence_time,
                        ref_pop_size  = 100000,
                        gwas_pop_size = 2000,
                        mutation_rate = 5e-2,
                        Ne_A          = 1000):
  
  #Now we want to define our population parameters
  Ne_B = Ne_A/2 #Define the population split size (in terms of Ne's). This is the REFERENCE PANEL POP
  Ne_C = Ne_A/2 #Define the population split size (in terms of Ne's). This is the GWAS PANEL POP.

  #Define our population configurations in msprime
  population_configurations = [
        msprime.PopulationConfiguration(
            #Source population (id = 0), goes extinct at time t
            sample_size  = 0,
            initial_size = Ne_A),
        msprime.PopulationConfiguration(
            sample_size  = int(ref_pop_size),
            initial_size = Ne_B), #Reference panel size
        msprime.PopulationConfiguration(
            sample_size  = int(gwas_pop_size),
            initial_size = Ne_C) #GWAS panel size
  ]
  #Define the demographic events
  demographic_events = [
    # Merging of pops B into A at some time t
    msprime.MassMigration(
        time         = divergence_time,
        source       = 2,
        destination  = 0,
        proportion   = 1.0),
    # Merging of pops C into A at some time t
    msprime.MassMigration(
        time         = divergence_time,
        source       = 1,
        destination  = 0,
        proportion   = 1.0),
    #No need to kill of Pop A, since we sample B & C; can keep it for future use
    # msprime.PopulationParametersChange(
    #     time=merge_time, initial_size=0, population_id=0)
  ]
  
  #Generate the recombination rate that we care about (here its no recomb/hotspot/no recomb)
  
  map_positions = [i*1 for i in range(0, 3)]
  map_rates = [0,1e-1,0]
  my_map = msprime.RecombinationMap(map_positions, map_rates)
  
  #Perform the actual msprime tree generation now
  result =          msprime.simulate(population_configurations = population_configurations,
                                     demographic_events        = demographic_events,
                                     mutation_rate             = mutation_rate,
                                     recombination_map         = my_map)
  return(np.asarray(allel.HaplotypeArray(result.genotype_matrix())))
                                     
                                     
                                     
            
