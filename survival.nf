/*

##############################################################################
##############################################################################

 TCGA survival pipeline:

 This pipeline allows us to calculate survivial curves related
 to individual gene segregation for individual cancers and the
 TCGA project overall.



##############################################################################
##############################################################################

*/

import java.text.SimpleDateFormat

// Define an output directory based on a timestamp:
def date = new java.util.Date()
def sdf = new SimpleDateFormat("yyMMdd.HHmmss")
def timestamp =  sdf.format(date)
output_dir = params.output_dir + "/" + timestamp + "/" + params.gene + "_" + params.low_q + "_" + params.high_q


process extract_count_matrices {

    tag "Extract count matrices"
    publishDir "${output_dir}/raw_counts", mode:"copy"
    container "ghcr.io/xyonetx/tcga-pipeline/pandas"
    cpus 2
    memory '7 GB'

    input:
        path(hdf5)
        path(ann)

    output:
        path("*.raw_counts.all.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/extract_all_count_matrices.py \
            -f ${hdf5} \
            -a ${ann} \
            -o raw_counts.all.tsv
        """
}


process segregate_by_expression {

    tag "Segregate expression on $raw_counts"
    publishDir "${output_dir}/annotations", mode:"copy", pattern: "*annotations*"
    publishDir "${output_dir}/normalized_counts", mode:"copy", pattern: "*.norm_counts.tsv"
    container "ghcr.io/xyonetx/tcga-pipeline/deseq2"
    cpus 4
    memory '14 GB'

    input:
        path(raw_counts)
        path(full_annotations)

    output:
        path("${tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv")
        path("${tcga_type}.norm_counts.tsv")

    script:
        tcga_type = raw_counts.name.split('\\.')[0]
        """
        /opt/software/miniconda/envs/deseq2/bin/Rscript \
            /opt/software/scripts/normalize_and_segregate.R \
            ${raw_counts} \
            ${full_annotations} \
            ${params.gene} \
            ${params.low_q} \
            ${params.high_q} \
            ${tcga_type}.annotations.${params.gene}_${params.low_q}_${params.high_q}.tsv \
            ${tcga_type}.norm_counts.tsv
        """
}

process calculate_individual_survival {
    tag "Calculate Kaplan-Meier estimates for individual type"
    publishDir "${output_dir}/kaplan_meier", mode:"copy"
    container "ghcr.io/xyonetx/tcga-pipeline/pandas"
    cpus 1
    memory '4 GB'

    input:
        path(ann)

    output:
        path('kaplan_meier.*.png')

    script:
        """
        /usr/bin/python3 /opt/software/scripts/calculate_survival.py ${ann}
        """
}

process calculate_overall_survival {
    tag "Calculate Kaplan-Meier estimates for entire cohort"
    publishDir "${output_dir}/kaplan_meier", mode:"copy"
    container "ghcr.io/xyonetx/tcga-pipeline/pandas"
    cpus 1
    memory '4 GB'

    input:
        path 'ann??.tsv'

    output:
        path('kaplan_meier.*.png')

    script:
        """
        /usr/bin/python3 /opt/software/scripts/calculate_survival.py ann*.tsv
        """
}

process plot_individual_count_distribution {
    tag "Plot count distribution for individual type"
    publishDir "${output_dir}/distribution_plots", mode:"copy"
    container "ghcr.io/xyonetx/tcga-pipeline/pandas"
    cpus 1
    memory '2 GB'

    input:
        path(nc)

    output:
        path('nc*.png')

    script:
        """
        /usr/bin/python3 /opt/software/scripts/plot_count_distribution.py \
            -i ${nc} \
            -g ${params.gene}
        """
}

workflow {
    raw_counts_ch = extract_count_matrices(params.hdf5, params.full_annotations).flatten()
    subtype_ann_ch, nc_ch = segregate_by_expression(raw_counts_ch, params.full_annotations)
    calculate_individual_survival(subtype_ann_ch)
    merged_ann_ch = subtype_ann_ch.collect()
    calculate_overall_survival(merged_ann_ch)
    plot_individual_count_distribution(nc_ch)
}