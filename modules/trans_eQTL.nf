/*
* STEP 6B - Run QTLtools in permutation mode
*/

process run_permutation_trans {
        tag "${qtl_subset} - ${batch_index}/${params.n_batches}"
        publishDir "${params.outdir}/trans/${qtl_subset}", mode: 'copy'
        container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'


        input:
        each batch_index
        tuple val(qtl_subset), file(bed), file(bed_index), file(fastqtl_bed), file(fastqtl_bed_index), file(vcf), file(vcf_index), file(covariate)

        output:
        tuple val(qtl_subset), file("${qtl_subset}.permutation.batch.${batch_index}.${params.n_batches}.best.txt.gz"), file("${qtl_subset}.permutation.batch.${batch_index}.${params.n_batches}.bins.txt.gz"), file("${qtl_subset}.permutation.batch.${batch_index}.${params.n_batches}.hits.txt.gz")

        script:
        """
        QTLtools trans --vcf $vcf --bed $bed --cov $covariate --chunk $batch_index ${params.n_batches} --out ${qtl_subset}.permutation.batch.${batch_index}.${params.n_batches} --window ${params.cis_window} --permute ${params.n_permutations}
        """

}

/*
 * STEP 7B - Merge permutation batches from QTLtools
 */
process merge_permutation_batches_trans {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/trans", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file('*.best.txt.gz'), file('*.bins.txt.gz'), file('*.hits.txt.gz')
    output:
    tuple val(qtl_subset), file("${qtl_subset}.all_best.gz"), file("${qtl_subset}.all_bins.gz"), file("${qtl_subset}.all_hits.gz")

    script:
    """
    ls *best.txt.gz | xargs cat > ${qtl_subset}.all_best.gz &
    ls *bins.txt.gz | xargs cat > ${qtl_subset}.all_bins.gz &
    ls *hits.txt.gz | xargs cat > ${qtl_subset}.all_hits.gz
    """
}

/*
 * STEP 8B - Run QTLtools in nominal mode
 */
process run_nominal_trans {
    tag "${qtl_subset} - ${batch_index}/${params.n_batches}"
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    each batch_index
    tuple val(qtl_subset), file(bed), file(bed_index), file(fastqtl_bed), file(fastqtl_bed_index), file(vcf), file(vcf_index), file(covariate)

    output:
    tuple val(qtl_subset), file("${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}.best.txt.gz"), file("${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}
.bins.txt.gz"), file("${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches}.hits.txt.gz")

    script:
    """
    QTLtools trans --vcf $vcf --bed $bed --cov $covariate --chunk $batch_index ${params.n_batches} --out ${qtl_subset}.nominal.batch.${batch_index}.${params.n_batches} --window ${
params.cis_window} --nominal
    """
}


/*
 * STEP 9B - Merge nominal batches from QTLtools
 */
process merge_nominal_batches_trans {
    tag "${qtl_subset}"
    publishDir "${params.outdir}/trans", mode: 'copy'
    container = 'quay.io/eqtlcatalogue/qtlmap:v20.05.1'

    input:
    tuple val(qtl_subset), file('*.best.txt.gz'), file('*.bins.txt.gz'), file('*.hits.txt.gz')
    output:
    tuple val(qtl_subset), file("${qtl_subset}.all_best.gz"), file("${qtl_subset}.all_bins.gz"), file("${qtl_subset}.all_hits.gz")

    script:
    """
    ls *best.txt.gz | xargs cat > ${qtl_subset}.all_best.gz &
    ls *bins.txt.gz | xargs cat > ${qtl_subset}.all_bins.gz &
    ls *hits.txt.gz | xargs cat > ${qtl_subset}.all_hits.gz
    """
}
