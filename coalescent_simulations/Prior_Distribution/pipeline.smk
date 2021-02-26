"""
#########README##########
Running the inference of coalescent simulations on the cluster
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations
We run this file on the directory:  /well/mcvean/mtutert/thesis/coalescent_coverage
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 100 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely, besides changing around default
"""
import itertools as itertools
print("Executing snakefile")
#Hard code in number of replicates for now...but make this flexible later on!
replicates = list(range(0,10)) #Current max of 25
replicates = [str(i) for i in replicates]
#Number of chunks (each file, total number of samples is chunks * nSamples)
chunks     = list(range(1,10))
chunks     = [str(i) for i in chunks]
#Divergence Times
divergence_times = ['0','10','50','125']
#Fst Values
fst_parameters   = ['0.0012','0.0052','0.0211','0.0506']
#Calculate all pairs of these values (note that the diagonal is what we want in the end)
paired_values    = list(itertools.product(divergence_times,fst_parameters))
#Convert this to a string that looks appropriate
def join_tuple_string(strings_tuple) -> str:
   return '_'.join(strings_tuple)

# joining all the tuples
result = map(join_tuple_string, paired_values)

# converting and printing the result
paired_values = (list(result))
print(paired_values) #First value is ithe divergence time, second value is the Fst used in the model

rule all:
    input:
        expand("results/pop_split/{paired_values}_split_LD_pooled_quantile_counts.RData", paired_values = paired_values),
        expand("results/pop_split/{paired_values}_split_AF_pooled_quantile_counts.RData", paired_values = paired_values),
rule draw_weights:
    output:
        "results/pop_split/{paired_values}_split_LD_Replicate_{replicates}_Chunk_{chunks}.RData",
        "results/pop_split/{paired_values}_split_AF_Replicate_{replicates}_Chunk_{chunks}.RData"
    params:
        replicates        = "{replicates}",
        chunk             = "{chunks}",
        paired_values     = "{paired_values}"
    script:
        "draw_weights.R"

rule calculate_quantiles:
    input:
        expand("results/pop_split/{paired_values}_split_LD_Replicate_{replicates}_Chunk_{chunks}.RData",replicates = replicates, chunks = chunks, paired_values = paired_values),
        expand("results/pop_split/{paired_values}_split_AF_Replicate_{replicates}_Chunk_{chunks}.RData",replicates = replicates, chunks = chunks, paired_values = paired_values)
    output:
        "results/pop_split/{paired_values}_split_LD_quantile_counts_replicate_{replicates}.RData",
        "results/pop_split/{paired_values}_split_AF_quantile_counts_replicate_{replicates}.RData"
    params:
        replicates        = "{replicates}",
        paired_values     = "{paired_values}"
    script:
        "calculate_quantiles.R"

rule pool_quantiles:
    input:
        expand("results/pop_split/{paired_values}_split_LD_quantile_counts_replicate_{replicates}.RData", replicates = replicates, paired_values = paired_values),
        expand("results/pop_split/{paired_values}_split_AF_quantile_counts_replicate_{replicates}.RData", replicates = replicates, paired_values = paired_values)
    output:
        "results/pop_split/{paired_values}_split_LD_pooled_quantile_counts.RData",
        "results/pop_split/{paired_values}_split_AF_pooled_quantile_counts.RData"
    params:
        paired_values     = "{paired_values}"
    script:
        "pool_quantiles.R"
# #
# rule plot_results:
#     input:
#         expand("results/pop_split/{paired_values}_split_LD_pooled_quantile_counts.RData", paired_values = paired_values),
#         expand("results/pop_split/{paired_values}_split_LD_pooled_quantile_counts.RData", paired_values = paired_values)
#     output:
#         "results/grid_qq_plot.jpeg
#     script:
#         "plot_results.R"
