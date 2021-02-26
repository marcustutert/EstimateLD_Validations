"""
#########README##########
Running the inference of OOA coalescent simulations on the cluster
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations
We run this file on the directory:  /well/mcvean/mtutert/thesis/coalescent_coverage
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 100 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely, besides changing around default
"""
###########################
#Define snakemakeparameters
###########################

import itertools as itertools
replicates = list(range(1,2)) #Current max of 25, figure how to bring this in directly
replicates = [str(i) for i in replicates]
#Number of chunks (each file, total number of samples is chunks * nSamples)
chunks     = list(range(1,2))
chunks     = [str(i) for i in chunks]
#Need to write out all possible combo's for the reference and gwas panel
#Three populations are YRI, CEU & CHB
populations = ["CEU","CHB","YRI"]
fst         = ["0.15","0.25","0.35"]
#Get all permutations of the GWAS and reference population_pairs1
#Will sort out the matching of them later on
population_pairs = []
for x in itertools.permutations(populations, 2):
    population_pairs.append('%s_GWAS_%s_Ref' %(x[0], x[1]))

#Calculate all pairs of these values (note that the diagonal is what we want in the end)
paired_values    = list(itertools.product(population_pairs,fst))
#Convert this to a string that looks appropriate
def join_tuple_string(strings_tuple) -> str:
   return '_'.join(strings_tuple)

# joining all the tuples
result = map(join_tuple_string, paired_values)

# converting and printing the result
paired_values = (list(result))
print(paired_values) #This has a format of: GWASPOP_GWAS_REFPOP_Ref_FST


##############################################
#Rules

rule all:
    input:
        expand("results/OOA/{paired_values}_AF_pooled_quantile_counts.RData", paired_values = paired_values),
        expand("results/OOA/{paired_values}_AF_pooled_quantile_counts.RData", paired_values = paired_values)


rule draw_weights:
    output:
        "results/OOA/{paired_values}_AF_Replicate_{replicates}_Chunk_{chunks}.RData",
        "results/OOA/{paired_values}_LD_Replicate_{replicates}_Chunk_{chunks}.RData"
    params:
        paired_values    = "{paired_values}",
        replicates       = "{replicates}",
        chunks           = "{chunks}"
    script:
        "OOA_draw_weights.R"

rule calculate_quantiles:
    input:
        expand("results/OOA/{paired_values}_LD_Replicate_{replicates}_Chunk_{chunks}.RData", replicates = replicates, chunks = chunks, paired_values = paired_values),
        expand("results/OOA/{paired_values}_AF_Replicate_{replicates}_Chunk_{chunks}.RData", replicates = replicates, chunks = chunks, paired_values = paired_values),
    output:
        "results/OOA/{paired_values}_LD_quantile_counts_replicate_{replicates}.RData",
        "results/OOA/{paired_values}_AF_quantile_counts_replicate_{replicates}.RData"
    params:
        replicates        = "{replicates}",
        paired_values     = "{paired_values}"
    script:
        "OOA_calculate_quantiles.R"

rule pool_quantiles:
    input:
        expand("results/OOA/{paired_values}_LD_quantile_counts_replicate_{replicates}.RData", replicates = replicates, paired_values = paired_values),
        expand("results/OOA/{paired_values}_AF_quantile_counts_replicate_{replicates}.RData", replicates = replicates, paired_values = paired_values)
    output:
        "results/OOA/{paired_values}_AF_pooled_quantile_counts.RData",
        "results/OOA/{paired_values}_LD_pooled_quantile_counts.RData"
    params:
        paired_values  = "{paired_values}"
    script:
        "OOA_pool_quantiles.R"

# rule plot_results:
#     input:
#         expand("results/OOA/{population_pairs}_AF_pooled_quantile_counts.RData", population_pairs = population_pairs),
#         expand("results/OOA/{population_pairs}_AF_pooled_quantile_counts.RData", population_pairs = population_pairs)
#     output:
#         "results/OOA/figures/CEU_GWAS_CHB_Ref_AF_result.jpeg"
#     script:
#         "OOA_plot_results.R"
