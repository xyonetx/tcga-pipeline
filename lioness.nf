params.exp_mtx = 's3://xyone-demo-data/exp_with_header.tsv'
params.motifs = 's3://xyone-demo-data/motifs.txt'
params.ppi = 's3://xyone-demo-data/ppi.txt'
params.ann = 's3://xyone-demo-data/demo_annotations.tsv'
params.output_dir = 's3://xyone-batch-results'

process run_lioness {
    tag "Run LIONESS expression on $exp_mtx"
    publishDir "${params.output_dir}/lioness", mode:"copy"
    container "docker.io/blawney/pandas"
    cpus 4
    memory '12 GB'

    input:
        path(exp_mtx)
        path(motifs)
        path(ppi)
        path(ann)

    output:
        path("lioness_output/*.csv")
        path("lioness_weights.genes.merged.tsv")
        path("lioness_weights.tf.merged.tsv")
    script:
        """
        /usr/bin/python3 /opt/software/scripts/run_lioness.py \
            -i ${exp_mtx} \
            -m ${motifs} \
            -p ${ppi} \
            -a ${ann} \
            -n ${task.cpus}

        /usr/bin/python3 /opt/software/scripts/merge_lioness.py \
           -d lioness_output \
           -g lioness_weights.genes.merged.tsv \
           -t lioness_weights.tf.merged.tsv
        """
}


workflow {
    exp_mtx_ch = Channel.fromPath(params.exp_mtx)
    motifs_ch = Channel.fromPath(params.motifs)
    ppi_ch = Channel.fromPath(params.ppi)
    ann_ch = Channel.fromPath(params.ann)
    run_lioness(exp_mtx_ch, motifs_ch, ppi_ch, ann_ch)
}