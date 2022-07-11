import os
import pandas as pd


# step:

##* genome signatures should exist in same folder as genomes.
#do: 
#- ANI of genome signature vs sgc unitigs sig (recovered using that genome as a query)

out_dir = config.get("output_dir", "output.sgc_ani")
logs_dir = os.path.join(out_dir, "logs")

sgc_dir = "/home/ctbrown/2020-hu-rerun"
SUBSETS = ["denticola", "gingivalis", "bacteriodes"]

# read genome lists
genomeD = {}
for subset in SUBSETS:
    genome_list = [x.strip() for x in open(f"{subset}_genomes.txt")]
    genomeD[subset] = genome_list

rule all:
    input: expand(f"{out_dir}/{{subset}}.ani.csv.gz", subset=SUBSETS)

rule sourmash_ANI:
    input: 
       ref_genome=f"{sgc_dir}/{{subset}}/{{genome}}.sig",
       sgc_unitigs=f"{sgc_dir}/podarV_k31_r1_search_oh0_{{subset}}/{{genome}}.contigs.sig",
    output:
       prefetch=f"{out_dir}/{{subset}}/{{genome}}.contigs.sig"
    log: f"{logs_dir}/{{subset}}/{{genome}}.prefetch.log"
    benchmark: f"{logs_dir}/{{subset}}/{{genome}}.prefetch.benchmark"
    conda: "conf/env/sourmash4.4.1.yml"
    shell:
        """
        sourmash prefetch {input.ref_genome} {input.sgc_unitigs} -o {output} 2> {log}
        """
        # could use compare, but prefetch gives more info  (contain, jaccard, ANI all at once)
# then plot?

rule aggregate_ANI_by_subset:
    input: lambda w: expand(f"{out_dir}/{{subset}}/{{genome}}.contigs.sig", subset = w.subset, genome=genomeD[w.subset])
    output: f"{out_dir}/{{subset}}.ani.csv.gz"
    run:
        # aggreate all prefetch csvs --> single csv
        aggDF = pd.concat([pd.read_csv(str(x), sep=",") for x in input])
        aggDF.to_csv(str(output), index=False)
            

