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
output_dir = params.output_dir + "/" + timestamp


def print_help() {
    log.info"""

    Hello world
    ${output_date}
    """.stripIndent()
}


process extract_tcga_type {

    tag "Extract tcga type: $tcga_type"
    publishDir "${output_dir}/${tcga_type}/raw_counts", mode:"copy"
    container "blawney/pandas"
    cpus 4
    memory '8 GB'

    input:
        path(hdf5)
        val tcga_type

    output:
        path("${tcga_type}.tsv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/extract_tcga_type.py \
            -f ${hdf5} \
            -t ${tcga_type} \
            -o ${tcga_type}.tsv
        """
}

process normalize_counts {

    tag "Normalize counts on $input_file"
    publishDir "${output_dir}/${tcga_type}/full_normalized_counts", mode:"copy"
    container "blawney/deseq2"
    cpus 4
    memory '8 GB'

    input:
        path(raw_counts)
        val tcga_type

    output:
        path("${tcga_type}.normalized.tsv")

    script:
        """
        /opt/software/miniconda/envs/deseq2/bin/Rscript /opt/software/normalize.R \
            ${raw_counts} \
            ${tcga_type}.normalized.tsv 
        """
}


workflow {

    if (params.help){
        print_help()
        exit 0
    }

    raw_count_ch= extract_tcga_type(params.hdf5, params.tcga_type)
    norm_count_ch = normalize_counts(raw_count_ch, params.tcga_type)

}