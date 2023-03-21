process run_lioness {
    tag "Run LIONESS expression on $exp_mtx"
    publishDir "${params.output_dir}/lioness", mode:"copy"
    container "blawney/pandas"
    cpus 4
    memory '8 GB'

    input:
        path(exp_mtx)
        path(motifs)
        path(ppi)

    output:
        path("lioness_output/*.csv")

    script:
        """
        /usr/bin/python3 /opt/software/scripts/run_lioness.py \
            -i ${exp_mtx} \
            -m ${motifs} \
            -p ${ppi} \
            -n ${task.cpus}
        """
}


workflow {
    exp_mtx_ch = Channel.fromPath('/Users/brianlawney/fpt/muc1_updated/toydata/exp_mtx_with_header.tsv')
    motifs_ch = Channel.fromPath('/Users/brianlawney/fpt/muc1_updated/toydata/motifs.txt')
    ppi_ch = Channel.fromPath('/Users/brianlawney/fpt/muc1_updated/toydata/ppi.txt')
    run_lioness(exp_mtx_ch, motifs_ch, ppi_ch)
}