/*

##############################################################################
##############################################################################

This process calculates edge weight differential targeting 
for an already-prepared LIONESS matrix.


##############################################################################
##############################################################################

*/

import java.text.SimpleDateFormat

// Define an output directory based on a timestamp:
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd.HHmmss")
def timestamp =  sdf.format(date)
output_dir = params.output_dir + "/" + timestamp + "/" + params.gene + "_" + params.low_q + "_" + params.high_q


process calc_target_scores {

    tag "Calculate target scores and stats"
    publishDir "${output_dir}/target_score_stats", mode:"copy"
    container "ghcr.io/xyonetx/tcga-pipeline/pandas"
    cpus 2
    memory '6 GB'

    input:
        path(lioness_file)
        path(ann)

    output:
        path("lioness_gene_targeting_scores.tsv")
        path("lioness_gene_targeting.mw_test.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/lioness_gene_target_scores.py \
            -i ${lioness_file} \
            -a ${ann} \
            -s lioness_gene_targeting_scores.tsv \
            -t lioness_gene_targeting.mw_test.tsv
        """
}



workflow {
    //lioness_file_ch = Channel.fromPath('s3://xyonetx-tcga-pipeline/results/ENSG00000185499_0.3_0.7/TCGA-PAAD/lioness_merged/TCGA-PAAD.ENSG00000185499_0.3_0.7.lioness_merged.csv')
    //ann_ch = Channel.fromPath('s3://xyonetx-tcga-pipeline/results/ENSG00000185499_0.3_0.7/TCGA-PAAD/annotations/TCGA-PAAD.annotations.ENSG00000185499_0.3_0.7.tsv')
    lioness_file_ch = Channel.fromPath('s3://xyone-demo-data/dummy/dummy_lioness.csv')
    ann_ch = Channel.fromPath('s3://xyone-demo-data/dummy/ann.tsv')
    calc_target_scores(lioness_file_ch, ann_ch)
}