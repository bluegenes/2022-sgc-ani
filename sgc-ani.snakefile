import os
import pandas as pd


# step:

##* genome signatures should exist in same folder as genomes.
#do: 
#- ANI of genome signature vs sgc unitigs sig (recovered using that genome as a query)

out_dir = config.get("output_dir", "output.sgc_ani")
logs_dir = os.path.join(out_dir, "logs")

sgc_dir = "/home/ctbrown/2020-rerun-hu"
SUBSETS = ["denticola", "gingivalis", "bacteroides"]

# read genome lists
genomeD = {}
for subset in SUBSETS:
    genome_list = [x.strip() for x in open(f"{subset}_genomes.txt")]
    genomeD[subset] = genome_list

rule all:
    input: expand(f"{out_dir}/{{subset}}.ani.csv", subset=SUBSETS)

rule sourmash_ANI:
    input: 
       ref_genome=f"{sgc_dir}/{{subset}}/{{genome}}.sig",
       sgc_unitigs=f"{sgc_dir}/podarV_k31_r1_search_oh0_{{subset}}/{{genome}}.contigs.sig",
    output:
        f"{out_dir}/{{subset}}/{{genome}}.prefetch.csv"
    log: f"{logs_dir}/{{subset}}/{{genome}}.prefetch.log"
    benchmark: f"{logs_dir}/{{subset}}/{{genome}}.prefetch.benchmark"
    conda: "conf/env/sourmash4.4.1.yml"
    shell:
        """
        sourmash prefetch --threshold-bp 0 {input.ref_genome} {input.sgc_unitigs} -o {output} 2> {log}
        """
        # could use compare, but prefetch gives more info  (contain, jaccard, ANI all at once)
# then plot?

rule write_prefetch_filelist:
    input: lambda w: expand(f"{out_dir}/{w.subset}/{{genome}}.prefetch.csv", genome=genomeD[w.subset])
    output: f"{out_dir}/{{subset}}.prefetch.filelist"
    run:
        with open(str(output), 'w') as outF:
            for inF in input:
                fn = os.path.abspath(str(inF))
                outF.write(f"{fn}\n")


rule aggregate_ANI_by_subset:
    input: f"{out_dir}/{{subset}}.prefetch.filelist",
    output: f"{out_dir}/{{subset}}.ani.csv"
    conda: "conf/env/sourmash4.4.1.yml"
    shell:
        """
        python recalc-ani.py --from-file {input} -o {output} 
        """
