"""
#########README##########
Running distributions of LD and associated code
Assume we have generated the coalescent simulations locally with files located in: /Users/marcustutert/Desktop/Oxford_Dphil/InferLD_Validations/coalescent_simulations
We run this file on the directory:  /well/mcvean/mtutert/thesis/distr_ld
With the command:
snakemake --snakefile PATH/TO/SNAKEMAKE-j 100 --max-status-checks-per-second 0.01 --profile /well/mcvean/mtutert/snakemake/profile -f --rerun-incomplete -n
Note that we may need to add specific things to the profile as we go on the fly--but unlikely, besides changing around default
"""
import itertools as itertools
print("Executing snakefile")
#Hard code in number of replicates for now...but make this flexible later on!
replicates = list(range(0,99)) #Current max of 25

rule all:
    input:
        #"msprime_data/gwas_sumstats"
        #"msprime_data/ref_98.csv"
        #expand("msprime_data/filtered_panels/ref_{replicates}", replicates = replicates)
        #expand("msprime_data/filtered_panels/filtered_ref_{replicates}", replicates = replicates)
        #expand("results/sumstat_imputed_ref_{replicates}", replicates = replicates)


#Generate msprime simulations
rule msprime_sims:
    output:
    #Last ref panel based on the replicates
        "msprime_data/ref_98.csv"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/python.yaml"
    script:
        "msprime_sims.py"

#Generate GWAS sumstats
rule generate_gwas_sumstats:
    output:
        "msprime_data/gwas_sumstats"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    script:
        "generate_gwas_sumstats.R"

#Generate list of filtered SNPs across reference panels
rule filtered_snps:
    output:
        "msprime_data/filtered_panels/ref_{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        replicates = "{replicates}"
    script:
        "filtered_snps.R"

#Match and extract snps
rule match_snps:
    output:
        "msprime_data/filtered_panels/filtered_ref_{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        replicates = "{replicates}"
    script:
        "match_snps.R"

#Run sumstat imputation
rule sumstat_imputation:
    output:
        "results/sumstat_imputed_ref_{replicates}"
    conda:
        "/well/mcvean/mtutert/snakemake/envs/r.yaml"
    params:
        replicates = "{replicates}"
    script:
        "sumstat_imputation.R"
