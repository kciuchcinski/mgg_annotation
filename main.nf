nextflow.enable.dsl=2

// ---------------------
// Initialise params
// ---------------------
params.PHAGES_DIR       = params.PHAGES_DIR       ?: "$baseDir/test/data"
params.OUTPUT_DIR       = params.OUTPUT_DIR       ?: "$baseDir/results"
params.INPUT_EXTENSION  = params.INPUT_EXTENSION  ?: "fasta"
params.PHAGE_MIN_LENGTH = params.PHAGE_MIN_LENGTH ?: 2000

params.BATCH_SIZE        = params.BATCH_SIZE        ?: 50
params.THREADS_PER_BATCH = params.THREADS_PER_BATCH ?: 8
params.SEARCH_TOOL       = params.SEARCH_TOOL       ?: "hhblits_omp"

// DB root
params.DB_ROOT = params.DB_ROOT ?: "$baseDir/databases"

// These three *must* be initialised before first use anywhere
params.GLIMMER_TRAIN = params.GLIMMER_TRAIN ?: "${params.DB_ROOT}/training-file_refseq.icm"

params.HHSUITE            = params.HHSUITE ?: [:]
params.HHSUITE.PHROGS     = params.HHSUITE.PHROGS ?: "${params.DB_ROOT}/phrogs_hhsuite_db/phrogs"
params.HHSUITE.PFAM       = params.HHSUITE.PFAM   ?: "${params.DB_ROOT}/pfamA_32/pfam"
params.HHSUITE.ECOD       = params.HHSUITE.ECOD   ?: "${params.DB_ROOT}/ECOD_F70/ECOD_F70_20200207"
params.HHSUITE.ALANDB     = params.HHSUITE.ALANDB ?: "${params.DB_ROOT}/AlanDavidson/profile-db/all_proteins"

params.METADATA                          = params.METADATA ?: [:]
params.METADATA.PHROGS_TABLE             = params.METADATA.PHROGS_TABLE     ?: "${params.DB_ROOT}/tables/phrog_annot_v4.tsv"
params.METADATA.ALAN_TABLE               = params.METADATA.ALAN_TABLE       ?: "${params.DB_ROOT}/tables/alan_annot.tsv"

// DB iterations – also need defaults
params.n_iter_phrogs = params.n_iter_phrogs ?: 1
params.n_iter_alan   = params.n_iter_alan   ?: 1
params.n_iter_pfam   = params.n_iter_pfam   ?: 2
params.n_iter_ecod   = params.n_iter_ecod   ?: 2

// Batch knobs
params.chunk_size_a3m = params.chunk_size_a3m ?: params.BATCH_SIZE
params.cpu_hhsuite    = params.cpu_hhsuite    ?: params.THREADS_PER_BATCH


Channel
  .value(file(params.PHAGES_DIR))
  .set { ch_phages_dir }

// 1) Preprocessing
process CHECK_INPUT {
    tag { "dir:${phage_dir.baseName}" }
    input:
    path phage_dir         // directory only
    output:
    path "phages.fasta",    emit: fasta
    path "metadata.tsv",    emit: metadata
    path "processed_input", emit: processed_dir
    shell:
    """
    set -euo pipefail
    mkdir -p processed_input

    # Call preprocessing with a directory, no file list interpolation.
    preprocessing \\
        --inputs-dir ${phage_dir} \\
        --ext ${params.INPUT_EXTENSION} \\
        --min-len ${params.PHAGE_MIN_LENGTH} \\
        --fasta phages.fasta \\
        --table metadata.tsv \\
        --outdir processed_input
    """
}

// 2) ORF calling
process PRODIGAL {
    input:
    path phages_fasta
    output:
    path "phages.prod"
    shell:
    """
    prodigal-gv -i ${phages_fasta} -o phages.prod -f sco -p meta -c
    """
}

process GLIMMER3 {
    input:
    path phages_fasta
    output:
    path "phages.predict"
    shell:
    """
    glimmer3 -o50 -g110 -t30 ${phages_fasta} ${params.GLIMMER_TRAIN} phages
    """
}

// 3) Concatenate and filter ORFs
process CONCAT_PRODIGAL {
  input:
  path prodigal_prod
  output:
  path "prodigal.csv"
  shell:
  """
  set -euo pipefail
  concat_prodigal --input ${prodigal_prod} --output prodigal.csv
  """
}

process CONCAT_GLIMMER {
  input:
  path glimmer_predict
  output:
  path "unfiltered-glimmer.csv"
  shell:
  """
  set -euo pipefail
  concat_glimmer --input ${glimmer_predict} --output unfiltered-glimmer.csv
  """
}

process FILTER_GLIMMER {
  input:
  path unfiltered_csv
  output:
  path "glimmer.csv", emit: glimmer_csv
  path "removed.csv"
  shell:
  """
  set -euo pipefail
  filter_glimmer --input ${unfiltered_csv} --output glimmer.csv --removed removed.csv
  """
}

// 4) Processing
process PROCESSING {
  input:
  path metadata_tsv
  path prodigal_csv
  path glimmer_csv
  output:
  path "ambigous.csv"
  path "orfs.csv"
  shell:
  """
  set -euo pipefail
  processing \
    --metadata ${metadata_tsv} \
    --prodigal ${prodigal_csv} \
    --glimmer ${glimmer_csv} \
    --ambigous ambigous.csv \
    --confident orfs.csv
  """
}

process EXTRACT_TRANSLATE_ORFS {
  input:
  path 'phages.fasta'
  path orfs_csv
  output:
  path "confident_orfs_tables/*.csv", emit: confident_tables
  path "erronous_orfs_tables/*.csv",  emit: erroneous_tables
  path "orfs/*.fasta",                emit: orfs_fasta
  path "proteins/*.fasta",            emit: proteins_fasta
  shell:
  """
  set -euo pipefail
  translation \
    --fasta phages.fasta \
    --orfs ${orfs_csv} \
    --outdir .
  """
}

process REORGANIZE_TABLE {
  input:
  path confident_orfs_tables
  output:
  path "confident_orfs.csv", emit: confident_orfs_csv
  shell:
  """
  set -euo pipefail
  reorganize_table --inputs *.csv --output confident_orfs.csv
  """
}

// 5) Concatenate proteins
process CONCAT_PROTEINS {
  input:
  path proteins_fasta
  output:
  path "proteins.fasta"
  shell:
  """
  set -euo pipefail
  concat_proteins *.fasta --output proteins.fasta
  """
}

// 6) Clustering
process CLUSTERING {
    input:
    path proteins_fasta
    output:
    path "raw_PCs.tsv", emit: raw_pcs
    path "raw_msa.a3m", emit: raw_msa
    shell:
    """
    mkdir -p tmp
    mmseqs createdb ${proteins_fasta} tmp/PROTEIN-DB > createdb.log 2>&1
    mmseqs cluster tmp/PROTEIN-DB tmp/CLUSTER-DB . --min-seq-id ${params.CLUSTERING.IDENTITY} -s ${params.CLUSTERING.SENSITIVITY} -c ${params.CLUSTERING.COVERAGE} -e ${params.CLUSTERING.EVAL} > cluster.log 2>&1
    mmseqs createtsv tmp/PROTEIN-DB tmp/PROTEIN-DB tmp/CLUSTER-DB raw_PCs.tsv > createtsv.log 2>&1
    mmseqs result2msa tmp/PROTEIN-DB tmp/PROTEIN-DB tmp/CLUSTER-DB tmp/CLU-MSA-DB --msa-format-mode 3 > result2msa.log 2>&1
    cp tmp/CLU-MSA-DB raw_msa.a3m
    """
}

process PCS2PROTEINS {
  input:
  path raw_pcs
  output:
  path "PCs2proteins.tsv"
  shell:
  """
  set -euo pipefail
  PCs2proteins --input ${raw_pcs} --output PCs2proteins.tsv
  """
}

process CLEAN_MSA {
  input:
  path raw_msa
  output:
  path "msa.a3m"
  shell:
  """
  set -euo pipefail
  clean_msa --input ${raw_msa} --output msa.a3m
  """
}

// 7) Split MSA
process SPLIT_MSA {
  input:
  tuple path(msa), path(pcs_map)
  output:
  path "msa_split/*.a3m", emit: per_pc_a3m
  shell:
  """
  set -euo pipefail
  split_msa --msa ${msa} --map ${pcs_map} --outdir msa_split
  """
}

// 8) Build FFindex
process BUILD_FFINDEX {
    tag "chunk_${task.hash.substring(0,8)}"
    
    input:
    // This tells Nextflow to put all files from the list into 'chunk/a3m/' automatically
    tuple val(pc_ids), path(a3m_files, stageAs: 'chunk/a3m/*') 

    output:
    tuple val(pc_ids), path("chunk/qdb.ffindex"), path("chunk/qdb.ffdata")

    script:
    """
    # Files are already in chunk/a3m/ because of stageAs
    cd chunk/a3m
    ffindex_build -s ../qdb.ffdata ../qdb.ffindex .
    """
}

// 9) Enrich via PHROGS
process ENRICH_MSA_PHROGS {
    tag "enrich_${task.hash.substring(0,8)}"
    cpus params.cpu_hhsuite
    input:
    tuple val(pc_ids), path(qidx), path(qdat)
    output:
    tuple val(pc_ids), path("enr/enr_a3m.ffindex"), path("enr/enr_a3m.ffdata")
    shell:
    """
    export OMP_STACKSIZE=32768
    ulimit -s unlimited || true

    mkdir -p enr
    ln -s \$(realpath ${qidx}) enr/qdb.ffindex
    ln -s \$(realpath ${qdat}) enr/qdb.ffdata
    cd enr
    ${params.SEARCH_TOOL} -i qdb -d ${params.HHSUITE.PHROGS} -oa3m enr_a3m -n 2 -cov 0 -p 0.95 -cpu ${task.cpus} > enrich.log 2>&1
    """
}

// 11) Batched DB searches (cartesian product db x enriched batches)
def dbMatrix = Channel.of(
    tuple("PHROGS", params.HHSUITE.PHROGS, params.n_iter_phrogs),
    tuple("ALANDB", params.HHSUITE.ALANDB, params.n_iter_alan),
    tuple("PFAM",   params.HHSUITE.PFAM,   params.n_iter_pfam),
    tuple("ECOD",   params.HHSUITE.ECOD,   params.n_iter_ecod)
)

process HHSUITE_SEARCH_BATCH {
    tag { "${dbname}_${task.hash.substring(0,8)}" }
    cpus params.cpu_hhsuite
    input:
    tuple val(dbname), val(dbpath), val(niter), val(pc_ids), path(ffidx), path(ffdat)
    output:
    tuple val(dbname), val(pc_ids), path("hhs/${dbname}.hhr.ffindex"), path("hhs/${dbname}.hhr.ffdata")
    shell:
    """
    export OMP_STACKSIZE=32768
    ulimit -s unlimited || true

    mkdir -p hhs
    ln -s \$(realpath ${ffidx}) hhs/enr_a3m.ffindex
    ln -s \$(realpath ${ffdat}) hhs/enr_a3m.ffdata
    cd hhs
    ${params.SEARCH_TOOL} -i enr_a3m -d ${dbpath} -o ${dbname}.hhr -cpu ${task.cpus} -mact 0.35 -p 50 -z 0 -v 0 -b 0 -qid 10 -cov 10 -E 1 -n ${niter} > ${dbname}.log 2>&1
    """
}

// 12) Unpack HHRs
process UNPACK_HHR {
    tag { "unpack_${dbname}_${task.hash.substring(0,8)}" }
    input:
    tuple val(dbname), val(pc_ids), path(ffidx), path(ffdat)
    output:
    tuple val(dbname), path(dbname)
    shell:
    """
    mkdir -p "${dbname}"
    ffindex_unpack ${ffdat} ${ffidx} ${dbname}
    """
}

process GATHER_HHR {
    tag "gather_${dbname}"

    input:
    tuple val(dbname), val(dirs)

    output:
    tuple val(dbname), path("${dbname}_hhr_files")

    shell:
    """
    mkdir -p "${dbname}_hhr_files"
    for d in ${dirs.join(' ')}; do
        if [ -d "\$d" ]; then
            find "\$d" -maxdepth 1 -type f -exec cp '{}' "${dbname}_hhr_files/" \\;
        fi
    done
    """
}

process COLLECT_HITS {
    if( params.WRITE_SEARCH_TABLE ) {
      publishDir "${params.OUTPUT_DIR}", mode: 'copy', pattern: 'search.tsv'
    }

    input:
    path pcs2proteins_tsv
    path phrogs_dir
    path alandb_dir
    path pfam_dir
    path ecod_dir
    path phrogs_tbl
    path alan_tbl

    output:
    path "search.tsv"

    script:
    """
    annotation_collect_hits \\
      --pcs2proteins ${pcs2proteins_tsv} \\
      --phrogs-dir ${phrogs_dir} \\
      --alan-dir   ${alandb_dir} \\
      --pfam-dir   ${pfam_dir} \\
      --ecod-dir   ${ecod_dir} \\
      --phrogs-table ${phrogs_tbl} \\
      --alan-table ${alan_tbl} \\
      --write-search ${params.WRITE_SEARCH_TABLE} \\
      --out-search search.tsv
    """
}

process FILTER_HITS {
    publishDir "${params.OUTPUT_DIR}", mode: 'copy', pattern: 'report.tsv'

    input:
    path search_tsv
    val phrogs_min_prob
    val phrogs_min_qcov
    val phrogs_min_tcov
    val domain_min_prob
    val domain_min_qcov
    val domain_min_tcov

    output:
    path "report.tsv"

    script:
    """
    annotation_filter_hits \\
      --search ${search_tsv} \\
      --phrogs-min-prob ${phrogs_min_prob} \\
      --phrogs-min-qcov ${phrogs_min_qcov} \\
      --phrogs-min-tcov ${phrogs_min_tcov} \\
      --domain-min-prob ${domain_min_prob} \\
      --domain-min-qcov ${domain_min_qcov} \\
      --domain-min-tcov ${domain_min_tcov} \\
      --out-report report.tsv
    """
}

process BUILD_ANNOTATION {
    if( params.WRITE_ANNOTATION_TABLE ) {
      publishDir "${params.OUTPUT_DIR}", mode: 'copy', pattern: 'annotation.tsv'
    }

    input:
    path confident_orfs_csv
    path pcs2proteins_tsv
    path report_tsv

    output:
    path "annotation_full.tsv", emit: annotation_full
    path "annotation.tsv", emit: annotation_public, optional: true

    script:
    """
    annotation_build_table \\
      --confident ${confident_orfs_csv} \\
      --pcs2proteins ${pcs2proteins_tsv} \\
      --report ${report_tsv} \\
      --write-annotation ${params.WRITE_ANNOTATION_TABLE} \\
      --out-annotation-full annotation_full.tsv \\
      --out-annotation-public annotation.tsv
    """
}

// 14) GenBank export
process GENBANK {
    publishDir "${params.OUTPUT_DIR}/genbanks", mode: 'copy', pattern: '*.gb'

    input:
    path annotation_full_tsv
    path metadata_csv

    output:
    path "*.gb"

    script:
    """
    genbank \\
      --annotation ${annotation_full_tsv} \\
      --metadata ${metadata_csv} \\
      --outdir .
    """
}

// --- Workflow Wiring ---
workflow {
    // Stage 1: Preprocessing & ORF Calling
    ch_check_input = CHECK_INPUT(ch_phages_dir)
    ch_prodigal    = PRODIGAL(ch_check_input.fasta)
    ch_glimmer     = GLIMMER3(ch_check_input.fasta)
    ch_concat_prod = CONCAT_PRODIGAL(ch_prodigal)
    ch_concat_glim = CONCAT_GLIMMER(ch_glimmer)

    ch_filter_glim = FILTER_GLIMMER(ch_concat_glim)
    ch_processing  = PROCESSING(ch_check_input.metadata, ch_concat_prod, ch_filter_glim.glimmer_csv)

    // Stage 2: Translation & Protein Prep
    ch_translate       = EXTRACT_TRANSLATE_ORFS(ch_check_input.fasta, ch_processing[1])
    ch_concat_proteins = CONCAT_PROTEINS(ch_translate.proteins_fasta)
    ch_reorganize      = REORGANIZE_TABLE(ch_translate.confident_tables)

    // Stage 3: Clustering & MSA
    ch_clustering   = CLUSTERING(ch_concat_proteins)
    ch_pcs2proteins = PCS2PROTEINS(ch_clustering.raw_pcs)
    ch_clean_msa    = CLEAN_MSA(ch_clustering.raw_msa)

    // Materialize paired inputs for SPLIT_MSA (no collect to avoid LinkedList payloads)
    def ch_msa_and_map = ch_clean_msa.combine(ch_pcs2proteins)
    ch_split_msa = SPLIT_MSA(ch_msa_and_map)

    // Stage 4: Batching & Enrichment
    // To avoid all big clusters (since they are passed in order), include a size-aware shuffle before chunking
    ch_a3m_chunks = ch_split_msa.per_pc_a3m
        .flatten()
        .map { f -> tuple(f.getBaseName(), f, f.size()) }   // add file size
        .toSortedList { a, b -> b[2] <=> a[2] }             // sort descending by size
        .flatMap()                                          // flatten back to items
        .map { it[0..1] }                                   // drop size, keep (id, file)
        .buffer(size: params.chunk_size_a3m, remainder: true)
        .map { chunk -> tuple(chunk.collect{ it[0] }, chunk.collect{ it[1] }) }


    ch_ffindex  = BUILD_FFINDEX(ch_a3m_chunks)
    ch_enriched = ENRICH_MSA_PHROGS(ch_ffindex)

    // Stage 5: DB Searches
    ch_search_batches = dbMatrix.combine(ch_enriched)
        .map { dbname, dbpath, niter, pc_ids, ffidx, ffdat ->
            tuple(dbname, dbpath, niter, pc_ids, ffidx, ffdat)
        }

    ch_hhsuite_batch = HHSUITE_SEARCH_BATCH(ch_search_batches)
    ch_unpacked_hhr  = UNPACK_HHR(ch_hhsuite_batch)

    ch_gathered_dirs = GATHER_HHR(ch_unpacked_hhr.groupTuple())

    // Create per-DB singleton channels for annotation from the gathered directories
    def ch_phrogs_dir = ch_gathered_dirs.filter { it[0] == 'PHROGS' }.map { it[1] }.collect()
    def ch_alandb_dir = ch_gathered_dirs.filter { it[0] == 'ALANDB' }.map { it[1] }.collect()
    def ch_pfam_dir   = ch_gathered_dirs.filter { it[0] == 'PFAM'   }.map { it[1] }.collect()
    def ch_ecod_dir   = ch_gathered_dirs.filter { it[0] == 'ECOD'   }.map { it[1] }.collect()

    def ch_phrogs_tbl = Channel.value(file(params.METADATA.PHROGS_TABLE))
    def ch_alan_tbl   = Channel.value(file(params.METADATA.ALAN_TABLE))

    // 6.1 Collect raw hits (search.tsv)
    ch_hits_raw = COLLECT_HITS(
        ch_pcs2proteins.collect(),
        ch_phrogs_dir,
        ch_alandb_dir,
        ch_pfam_dir,
        ch_ecod_dir,
        ch_phrogs_tbl,
        ch_alan_tbl,
    )

    // 6.2 Filter hits -> report.tsv
    ch_hits_filtered = FILTER_HITS(
        ch_hits_raw,
        Channel.value(params.ANNOTATION_FILTERS.PHROGS.MIN_PROB),
        Channel.value(params.ANNOTATION_FILTERS.PHROGS.MIN_QCOV),
        Channel.value(params.ANNOTATION_FILTERS.PHROGS.MIN_TCOV),
        Channel.value(params.ANNOTATION_FILTERS.DOMAIN_DBS.MIN_PROB),
        Channel.value(params.ANNOTATION_FILTERS.DOMAIN_DBS.MIN_QCOV),
        Channel.value(params.ANNOTATION_FILTERS.DOMAIN_DBS.MIN_TCOV)
    )

    // 6.3 Build per-ORF annotation table
    ch_annotation_full = BUILD_ANNOTATION(
        ch_reorganize.confident_orfs_csv,
        ch_pcs2proteins.collect(),
        ch_hits_filtered
    )

    // If ch_annotation_full captures the process result object:
    GENBANK(BUILD_ANNOTATION.out.annotation_full, ch_check_input.metadata)
}
