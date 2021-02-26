"""
#########README##########
Running the inference of coalescent simulations on the cluster
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations
We run this file on the directory:  /well/mcvean/mtutert/thesis/coalescent_coverage
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 100 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely
"""

print("Executing snakefile")
#Check how many replicates we have
replicates = list(range(0,10))
replicates = [str(i) for i in replicates]
#Check how many totaly samples we have
Sample_Index   = list(range(1,10)) #### DONT HARD CODE THIS IN
Sample_Index   = [str(i) for i in Sample_Index]

#
rule all:
    input:
        expand("inference_results/ld_results/ld_rep_{replicates}_sample_{Sample_Index}", replicates = replicates, Sample_Index = Sample_Index)

# rule grab_data:
#     input:
#         expand("inference_results/true_quantiles_rep_{replicates}",replicates = replicates)
#     output:
#         "inference_results/pooled_quantiles"
#     conda:
#         "/well/mcvean/mtutert/snakemake/envs/r.yaml"
#     script:
#         "pool_results.R"

rule install_package:
        script:
            "install_package.R"

rule coalescent_sims:
    output:
        "inference_results/inference_results_rep_{replicates}"
    params:
        prefix="{replicates}"
    script:
        "ld_inference.R"

rule calculate_ld_per_sample:
    input:
        "inference_results/inference_results_rep_{replicates}"
    output:
        "inference_results/ld_results/ld_rep_{replicates}_sample_{Sample_Index}"
    params:
        prefix = "{Sample_Index}",
        suffix = "{replicates}"
    script:
        "per_sample_ld_calculation.R"
#
# rule summarize_ld
#     input:
#         "inference_results/ld_results/ld_rep_{replicates}_sample_{Sample_Index}"
#     params:
#         prefix = "{Sample_Index}",
#         suffix = "{replicates}"
#     script:
#         "per_sample_ld_calculation.R"
#
#
