import numpy as np
import allel
import math as math

def marcus():
  return(print("hello"))

#Create a function that will generate 1000x replicates of msprime seeded trees, with split population (A->(B&C))) with
def ancestral_pop_split(divergence_time,
                        seed,
                        ref_pop_size,
                        gwas_pop_size,
                        mutation_rate,
                        length,
                        Ne_A, 
                        num_replicates,
                        statistics):

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


  #Perform the actual msprime tree generation now
  replicates = list(msprime.simulate(population_configurations = population_configurations,
                                     demographic_events        = demographic_events,
                                     mutation_rate             = mutation_rate,
                                     length                    = length,
                                     random_seed               = seed,
                                     num_replicates            = num_replicates,
                                     recombination_rate        = 1e-8) )

  if(statistics == 0): #Return replicates tree to perform Fst inference on (optional, Yan; ignore)
    return(replicates)

  if(statistics == 1): #Go through and return the actual haplotype array from GWAS & Ref population
    #Iterate through the trees
    trees = []
    for i, tree_sequence in enumerate(replicates):
      #Convert to haplotype array using scikit allel
      trees.append(np.asarray(allel.HaplotypeArray(tree_sequence.genotype_matrix())))

    b=[]
    for i in range(num_replicates):
        total_population = trees[i]
        total_population = total_population.transpose()
        #Split populations into two
        ref_pop          = total_population[0:ref_pop_size,]
        GWAS_pop         = total_population[ref_pop_size:,]
        b.append(GWAS_pop) #First list element
        b.append(ref_pop)  #Second list element
    #Remove list comprehension
    b = [b[i:i+2] for i in range(0, len(b), 2)]
    return(b)
