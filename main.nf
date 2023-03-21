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

// remove genes with cohort-wise means lower than this threshold:
// (where cohort is the subset of samples AFTER subsetting for our
// high and low expressed gene of interest)
params.min_reads = 5

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


process run_gsea {

    tag "Run GSEA on tcga type: $params.tcga_type"
    publishDir "${output_dir}/${params.tcga_type}/gsea", mode:"copy"
    container "blawney/gsea"
    cpus 4
    memory '8 GB'

    input:
        path(norm_counts)
        path(annotations)

    output:
        path("${gct_file}")
        path("${cls_file}")
        path("${params.tcga_type}.gsea_results.zip")

    script:
        def gct_file_template = "%s.%s_%s_%s.high_vs_low.gct"
        def cls_file_template = "%s.%s_%s_%s.high_vs_low.cls"
        gct_file = String.format(gct_file_template, params.tcga_type, params.gene, params.low_q, params.high_q)
        cls_file = String.format(cls_file_template, params.tcga_type, params.gene, params.low_q, params.high_q)
        """
        /usr/bin/python3 /opt/software/scripts/prep_files.py \
            -f ${norm_counts} \
            -a ${annotations} \
            -g ${gct_file} \
            -c ${cls_file} \
            -t ${params.min_reads}

        /opt/software/gsea/GSEA_4.3.2/gsea-cli.sh GSEA \
            -res "${gct_file}" \
            -cls "${cls_file}#high_versus_low" \
            -gmx ftp.broadinstitute.org://pub/gsea/msigdb/human/gene_sets/h.all.v2023.1.Hs.symbols.gmt \
            -chip ftp.broadinstitute.org://pub/gsea/msigdb/human/annotations/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -out /gsea/ \
            -rpt_label "${params.tcga_type}" \
            -zip_report true \
            -collapse Collapse \
            -mode Max_probe \
            -norm meandiv \
            -nperm 1000 \
            -permute phenotype \
            -rnd_seed timestamp \
            -rnd_type no_balance \
            -scoring_scheme weighted \
            -metric Signal2Noise \
            -sort real \
            -order descending \
            -create_gcts false \
            -create_svgs false \
            -include_only_symbols true \
            -make_sets true \
            -median false \
            -num 100 \
            -plot_top_x 20 \
            -save_rnd_lists false \
            -set_max 500 \
            -set_min 15

        /usr/bin/python3 /opt/software/scripts/move_final_files.py \
            -p "/gsea/${params.tcga_type}*/*.zip" \
            -o ${params.tcga_type}.gsea_results.zip
        """
}


process map_ensg_to_symbol {
    tag "Run ENSG to symbol gene mapping $exp_mtx"
    publishDir "${output_dir}/${params.tcga_type}/normalized_counts", mode:"copy"
    container "docker.io/blawney/pandas"
    cpus 2
    memory '4 GB'

    input:
        path(exp_mtx)

    output:
        path("${params.tcga_type}.deseq2_norm_counts.symbol_remapped.all.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/map_ensg_to_symbol.py \
            -i ${exp_mtx} \
            -m /opt/software/resources/Human_Ensembl_Gene_ID_MSigDB.v2023.1.Hs.chip \
            -o ${params.tcga_type}.deseq2_norm_counts.symbol_remapped.all.tsv
        """
}

workflow {

    if (params.help){
        print_help()
        exit 0
    }

    raw_count_ch= extract_tcga_type(params.hdf5)
    (dge_results_ch, norm_counts_ch, ann_ch) = run_dge(raw_count_ch)
    run_gsea(norm_counts_ch, ann_ch)
    norm_counts_symbol_remapped_ch = map_ensg_to_symbol(norm_counts_ch)

}