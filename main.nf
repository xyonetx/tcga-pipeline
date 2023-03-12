/*

##############################################################################
##############################################################################

 TCGA exploration pipeline:

 This pipeline allows us to parallelize a process over the entire TCGA
 project.



##############################################################################
##############################################################################

*/

import java.text.SimpleDateFormat

/*
    START default params:
*/

// unless specified, skip help message
params.help = false

/*
    END default params:
*/


// Define an output directory based on a timestamp:
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd.HHmmss")
def timestamp =  sdf.format(date)
output_dir = params.output_dir + "/" + timestamp + "/" + params.gene + "_" + params.low_q + "_" + params.high_q


def print_help() {
    log.info"""

    Hello world
    ${output_date}
    """.stripIndent()
}


process extract_tcga_type {

    tag "Extract tcga type: $params.tcga_type"
    publishDir "${output_dir}/${params.tcga_type}/raw_counts", mode:"copy"
    container "blawney/pandas"
    cpus 4
    memory '8 GB'

    input:
        path(hdf5)

    output:
        path("${params.tcga_type}.raw_counts.all.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/extract_tcga_type.py \
            -f ${hdf5} \
            -t ${params.tcga_type} \
            -o ${params.tcga_type}.raw_counts.all.tsv
        """
}

process run_dge {

    tag "Run differential expression on $raw_counts"
    publishDir "${output_dir}/${params.tcga_type}/normalized_counts", mode:"copy", pattern: "*.deseq2_norm_counts.all.tsv"
    publishDir "${output_dir}/${params.tcga_type}/dge_results", mode:"copy", pattern: "*.deseq2_results*"
    publishDir "${output_dir}/${params.tcga_type}/annotations", mode:"copy", pattern: "*.annotations*"
    container "blawney/deseq2"
    cpus 4
    memory '8 GB'

    input:
        path(raw_counts)

    output:
        path("${params.tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv")
        path("${params.tcga_type}.deseq2_norm_counts.all.tsv")
        path("${params.tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv")

    script:
        """
        /opt/software/miniconda/envs/deseq2/bin/Rscript /opt/software/scripts/deseq2.R \
            ${raw_counts} \
            ${params.gene} \
            ${params.low_q} \
            ${params.high_q} \
            ${params.tcga_type}.deseq2_results.${params.gene}_${params.low_q}_${params.high_q}.high_vs_low.tsv \
            ${params.tcga_type}.deseq2_norm_counts.all.tsv \
            ${params.tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv
        """
}


workflow {

    if (params.help){
        print_help()
        exit 0
    }

    raw_count_ch= extract_tcga_type(params.hdf5)
    (dge_results_ch, norm_counts_ch, ann_ch) = run_dge(raw_count_ch)

}