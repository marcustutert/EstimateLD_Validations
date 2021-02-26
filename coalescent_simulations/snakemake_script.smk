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
replicates = list(range(0,10))
replicates = [str(i) for i in replicates]
print(replicates)

rule all:
    input:
        "inference_results/pooled_quantiles"

rule grab_data:
    input:
        expand("inference_results/true_quantiles_rep_{replicates}",replicates = replicates)
    output:
        "inference_results/pooled_quantiles"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "pool_results.R"

rule coalescent_sims:
    output:
        "inference_results/true_quantiles_rep_{replicates}"
    params:
        prefix="{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "coverage_inference.R"
