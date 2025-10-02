"""
Phage oriented ORF prediction tool
by Janusz Koszucki, Jade Leconte, Wanangwa Ndovie, Rafal Mostowy

Phage-oriented pipeline for high-confidence Open Reading Frames prediction (genes)
Algorithm copied from Multiphage-2 tool.

version: 0.2
date: 20.07.2023
authors: Janusz Koszucki, Wanangwa Ndovie, Rafal Mostowy
"""

# Run with Apptainer
container: "phage_annotation.sif"

# load modules
from dependencies.scripts.utils import bcolors, checkpoint, get_n_iterations
from pathlib import Path
import pandas as pd
from Bio import SeqIO


# start
print(f"{bcolors.OKGREEN}Running Another Annotation Phage Tool! {bcolors.ENDC}", end='\n')

# warning
print(f"{bcolors.WARNING}WARNINGS!{bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}1. No solid checkpoints on input. {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}2. Optional flags to implement. Rename files, folders and variables when ugly. {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}3. Clean unused output {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}4. Works only when running all databases {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}5. Removing conda environments in some rules (HHsuite) can significantly speed up calculations {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}6. Remove redundant files (eg, 1_MSA_VIEW & 2_MSA_SEARCH) {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}7. SetupTools for setting up HHsuite databases automatically {bcolors.ENDC}", end='\n')
print(f"{bcolors.WARNING}8. ALANDB HHsuite database on DropBox is (?) corrupted (multiple HMM profiles are empty, due to inconsistent number of columns in alignments from Alan) {bcolors.ENDC}", end='\n')



# select databases
SEARCH_TOOL = config['SEARCH_TOOL']
DBs = ['PHROGS', 'ALANDB', 'PFAM', 'ECOD']


# paths & params
PHAGES_DIR = config['PHAGES_DIR']
OUTPUT_DIR = config['OUTPUT_DIR']
EXTENSION = config['INPUT_EXTENSION']
PHAGE_MIN_LENGTH = config['PHAGE_MIN_LENGTH']

IDENTITY = config['CLUSTERING']['IDENTITY']
COVERAGE = config['CLUSTERING']['COVERAGE']
EVAL = config['CLUSTERING']['EVAL']
SENSITIVITY = config['CLUSTERING']['SENSITIVITY']

GLIMMER_TRAIN = config['GLIMMER_TRAIN']

# output
IN_DIR_PROCESSED = Path(OUTPUT_DIR, '1_PROCESSED_INPUT')
ORF_PREDICTION_DIR = Path(OUTPUT_DIR, '2_ORF_PREDICTION')
ANNOTATION_DIR = Path(OUTPUT_DIR, '3_ANNOTATION')
GENBANK_DIR = Path(OUTPUT_DIR, '4_GENBANK')


### intermediate folders

# ORF prediction
PRODIGAL_DIR = Path(ORF_PREDICTION_DIR, '1_GENE_CALLING', '1_PRODIGALGV')
GLIMMER_DIR = Path(ORF_PREDICTION_DIR, '1_GENE_CALLING', '2_GLIMMER')

ORF_PROCESSING_DIR = Path(ORF_PREDICTION_DIR, '2_PROCESSING')
ORFS_DIR = Path(ORF_PREDICTION_DIR, '3_ORFS')
PROTEINS_DIR = Path(ORF_PREDICTION_DIR, '4_PROTEINS')


# annotation
VERSION=f'IDENT{str(int(IDENTITY * 100))}_COV{str(int(COVERAGE * 100))}'
CLUSTERING_DIR = Path(ANNOTATION_DIR, f'1_CLUSTERING_{VERSION}')
MSA_DIR = Path(ANNOTATION_DIR, '2_MSA')
HHSUITE_DIR = Path(ANNOTATION_DIR, '3_HHSUITE')

PHROGS = Path(config['HHSUITE']['PHROGS'])
PFAM = Path(config['HHSUITE']['PFAM'])
ECOD = Path(config['HHSUITE']['ECOD'])
ALANDB = Path(config['HHSUITE']['ALANDB'])

PHROGs_N_ITERATIONS, ALANDB_N_ITERATIONS, PFAM_N_ITERATIONS, ECOD_N_ITERATIONS = get_n_iterations(SEARCH_TOOL)

PHROGS_TABLE = Path(config['METADATA']['PHROGS_TABLE'])
MGG_PHROGS_TABLE = Path(config['METADATA']['MGG_PHROGS_TABLE'])
ALAN_TABLE = Path(config['METADATA']['ALAN_TABLE'])


                        #############################################
                        ######## PREPROCESSING & CHECKPOINTS ########
                        #############################################



# process input
discarded_phages = checkpoint(PHAGES_DIR, IN_DIR_PROCESSED, PHAGE_MIN_LENGTH, EXTENSION)
phages, = glob_wildcards(Path(IN_DIR_PROCESSED, '{phages}.' + EXTENSION))

# checkpoint
def trigger_search(wildcards):
    PCs_checkpoint = checkpoints.split_msa.get(**wildcards).output[0]
    global PCs
    PCs, = glob_wildcards(Path(PCs_checkpoint, '{PC}.a3m'))
    return expand(Path(HHSUITE_DIR, '{DB}', '{PC}.hhr'), DB=DBs, PC=PCs)



                    ################################
                    ############ TARGET ############
                    ################################

rule target:
    input:
        Path(ORF_PREDICTION_DIR, 'metadata.tsv'),               # preprocessed
        Path(PRODIGAL_DIR, 'prodigal.csv'),                     # run prodigal-gv
        Path(GLIMMER_DIR, 'glimmer.csv'),                       # run glimmer
        Path(ORF_PROCESSING_DIR, 'orfs.csv'),                   # processing results
        Path(ORF_PREDICTION_DIR, 'confident_orfs.csv'),         # orfs/proteins
        Path(ANNOTATION_DIR, 'proteins.fasta'),                 # concat proteins
        Path(ANNOTATION_DIR, 'msa.a3m'),
        Path(CLUSTERING_DIR, 'raw_PCs.tsv'),                    # expicitly ask for clustering
        Path(ANNOTATION_DIR, 'PCs2proteins.tsv'),               # map PCs to proteins
        Path(OUTPUT_DIR, 'annotation.tsv'),                     # annotation table
        Path(GENBANK_DIR)                                       # genbank per phage
        


        #####################################
        ######### CHECK INPUT FILES #########
        #####################################

# verify input
rule check_input:
    input: expand(Path(IN_DIR_PROCESSED, '{phage}.' + EXTENSION), phage=phages)
    output:
        fasta=Path(ORF_PREDICTION_DIR, 'phages.fasta'),
        table=Path(ORF_PREDICTION_DIR, 'metadata.tsv')
    script: 'dependencies/scripts/preprocessing.py'


        ################################
        ######## ORF PREDICTION ########
        ################################


### orf calling
# run prodigal (runs on single core)
rule prodigal:
    input: Path(ORF_PREDICTION_DIR, 'phages.fasta'),
    output: Path(PRODIGAL_DIR, 'phages.prod')
    threads: workflow.cores * 1.0
    shell: 'prodigal-gv -i {input} -o {output} -f sco -p meta -c'

# run glimmer (runs all cores)
rule glimmer3:
    input: 
        fasta=Path(ORF_PREDICTION_DIR, 'phages.fasta'),
        train=GLIMMER_TRAIN
    output: Path(GLIMMER_DIR, 'phages.predict'),
    threads: workflow.cores * 1.0
    shell:
        'glimmer3 -o50 -g110 -t30 {input.fasta} {input.train} {output}; '
        'mv {output}.predict {output}; '


### concatenate results from orf calling
# concat prodigal
rule concat_prodigal:
    input: Path(PRODIGAL_DIR, 'phages.prod')
    output: Path(PRODIGAL_DIR, 'prodigal.csv')
    script: 'dependencies/scripts/concat_prodigal.py'

# concat glimmer
rule concat_glimmer:
    input: Path(GLIMMER_DIR, 'phages.predict')
    output: Path(GLIMMER_DIR, 'unfiltered-glimmer.csv')
    script: 'dependencies/scripts/concat_glimmer.py'

# filter glimmer obvious false positives
rule filter_glimmer:
    input: Path(GLIMMER_DIR, 'unfiltered-glimmer.csv')
    output:
        Path(GLIMMER_DIR, 'glimmer.csv'),
        Path(GLIMMER_DIR, 'removed.csv')
    script:'dependencies/scripts/filter_glimmer.py'


### process results
# get high-confidence orfs
rule processing:
    input:
        Path(ORF_PREDICTION_DIR, 'metadata.tsv'),
        Path(PRODIGAL_DIR, 'prodigal.csv'),
        Path(GLIMMER_DIR, 'glimmer.csv')
    output:
        Path(ORF_PROCESSING_DIR, 'ambigous.csv'),
        Path(ORF_PROCESSING_DIR, 'orfs.csv')
    script: 'dependencies/scripts/processing.py'

# extract ORFs & translate to proteins
rule extract_translate_orfs:
    input:
        Path(IN_DIR_PROCESSED, '{phage}.' + EXTENSION),
        Path(ORF_PROCESSING_DIR, 'orfs.csv')
    output:
        Path(ORF_PROCESSING_DIR, 'confident_orfs_tables', '{phage}.csv'),
        Path(ORF_PROCESSING_DIR, 'erronous_orfs_tables', '{phage}.csv'),
        Path(ORFS_DIR, '{phage}.fasta'),
        Path(PROTEINS_DIR, '{phage}.fasta')
    script:'dependencies/scripts/translation.py'

# reorganize results in one table
rule reorganize_table:
    input: expand(Path(ORF_PROCESSING_DIR, 'confident_orfs_tables', '{phage}.csv'), phage=phages)
    output: Path(ORF_PREDICTION_DIR, 'confident_orfs.csv')
    script:'dependencies/scripts/reorganize_table.py'



        #######################################
        ######## FUNCTIONAL ANNOTATION ########
        #######################################

# concatenate proteins
rule concat_proteins:
    input: expand(Path(PROTEINS_DIR, '{phage}.fasta'), phage=phages)
    output: Path(ANNOTATION_DIR, 'proteins.fasta')
    script: 'dependencies/scripts/concat_proteins.py'


# cluster proteins
rule clustering:
    input: Path(ANNOTATION_DIR, 'proteins.fasta')
    output:
        mmseqs_dir=directory(Path(CLUSTERING_DIR)),
        clusters=Path(CLUSTERING_DIR, 'raw_PCs.tsv'),
        msa=Path(CLUSTERING_DIR, 'raw_msa.a3m'),
        protein_db=Path(CLUSTERING_DIR, 'tmp', 'PROTEIN-DB'),
        cluster_db=Path(CLUSTERING_DIR, 'tmp', 'CLUSTER-DB'),
        clust2msa=Path(CLUSTERING_DIR, 'tmp', 'CLU-MSA-DB')
    params:
        IDENTITY=IDENTITY,
        COVERAGE=COVERAGE,
        EVAL=EVAL,
        SENSITIVITY=SENSITIVITY
    shell:
        "mmseqs createdb {input} {output.protein_db} > {output.mmseqs_dir}/createdb.log 2>&1; "
        "mmseqs cluster {output.protein_db} {output.cluster_db} {output.mmseqs_dir} --min-seq-id {params.IDENTITY} -s {params.SENSITIVITY} -c {params.COVERAGE} -e {params.EVAL} > {output.mmseqs_dir}/cluster.log 2>&1; "
        "mmseqs createtsv {output.protein_db} {output.protein_db} {output.cluster_db} {output.clusters} > {output.mmseqs_dir}/createtsv.log 2>&1; "
        "mmseqs result2msa {output.protein_db} {output.protein_db} {output.cluster_db} {output.clust2msa} --msa-format-mode 3 > {output.mmseqs_dir}/result2msa.log 2>&1; "
        "cp {output.clust2msa} {output.msa}; " # copy results
        "touch {output.protein_db} {output.cluster_db} {output.clust2msa}; " # create snakemake dummy files exist


# map clusters to proteins
rule PCs2proteins:
    input: Path(CLUSTERING_DIR, 'raw_PCs.tsv')
    output: Path(ANNOTATION_DIR, 'PCs2proteins.tsv')
    script: 'dependencies/scripts/PCs2proteins.py'


# clean msa
rule clean_msa:
    input: Path(CLUSTERING_DIR, 'raw_msa.a3m')
    output: Path(ANNOTATION_DIR, 'msa.a3m')
    script: 'dependencies/scripts/clean_msa.py'


# msa into seperate files
### danger: using too much RAM 
### solution: increase threads to run less rules in parallel
checkpoint split_msa:
    input: 
        msa=Path(ANNOTATION_DIR, 'msa.a3m'),
        PCs2proteins=Path(ANNOTATION_DIR, 'PCs2proteins.tsv')
    output:
        msa4view=directory(Path(MSA_DIR, '1_MSA_VIEW')),
        msa4search=directory(Path(MSA_DIR, '2_MSA_SEARCH'))
    script: 'dependencies/scripts/split_msa.py'


# enrich msa with PHROGS
rule enrich_msa:
    input: Path(MSA_DIR, '2_MSA_SEARCH', '{PC}.a3m')
    output: Path(MSA_DIR, '3_MSA_PHROGS', '{PC}.a3m')
    params: PHROGS=PHROGS
    shell: 'hhblits -i "{input}" -d "{params.PHROGS}" -oa3m "{output}" -n 2 -cov 0 -p 0.95 > {output}.log 2>&1'


                #######################
                ######## SEARCH #######
                #######################

# search PHROGs (one iteration)
rule PHROGs:
    input: Path(MSA_DIR, '3_MSA_PHROGS', '{PC}.a3m')
    output: Path(HHSUITE_DIR, 'PHROGS', '{PC}.hhr')
    params: 
        SEARCH_TOOL=SEARCH_TOOL,
        PHROGS=PHROGS,
        PHROGs_N_ITERATIONS=PHROGs_N_ITERATIONS
    shell: '{params.SEARCH_TOOL} -i "{input}" -d "{params.PHROGS}" -o "{output}" -cpu 2 -mact 0.35 -p 50 -z 0 -v 0 -b 0 -qid 10 -cov 10 -E 1 {params.PHROGs_N_ITERATIONS} > {output}.log 2>&1'


# search ALANDB (one iteration)
rule ALANDB:
    input: Path(MSA_DIR, '3_MSA_PHROGS', '{PC}.a3m')
    output: Path(HHSUITE_DIR, 'ALANDB', '{PC}.hhr')
    params: 
        SEARCH_TOOL=SEARCH_TOOL,
        ALANDB=ALANDB,
        ALANDB_N_ITERATIONS=ALANDB_N_ITERATIONS
    shell: '{params.SEARCH_TOOL} -i "{input}" -d "{params.ALANDB}" -o "{output}" -cpu 2 -mact 0.35 -p 50 -z 0 -v 0 -b 0 -qid 10 -cov 10 -E 1 {params.ALANDB_N_ITERATIONS} > {output}.log 2>&1'


# search PFAM (two interations)
rule PFAM:
    input: Path(MSA_DIR, '3_MSA_PHROGS', '{PC}.a3m')
    output: Path(HHSUITE_DIR, 'PFAM', '{PC}.hhr')
    params: 
        SEARCH_TOOL=SEARCH_TOOL,
        PFAM=PFAM,
        PFAM_N_ITERATIONS=PFAM_N_ITERATIONS
    shell: '{params.SEARCH_TOOL} -i "{input}" -d "{params.PFAM}" -o "{output}" -cpu 2 -mact 0.35 -p 50 -z 0 -v 0 -b 0 -qid 10 -cov 10 -E 1 {params.PFAM_N_ITERATIONS} > {output}.log 2>&1'


# search ECOD (two interations)
rule ECOD:
    input: Path(MSA_DIR, '3_MSA_PHROGS', '{PC}.a3m')
    output: Path(HHSUITE_DIR, 'ECOD', '{PC}.hhr')
    params: 
        SEARCH_TOOL=SEARCH_TOOL,
        ECOD=ECOD,
        ECOD_N_ITERATIONS=ECOD_N_ITERATIONS
    shell: '{params.SEARCH_TOOL} -i "{input}" -d "{params.ECOD}" -o "{output}" -cpu 2 -mact 0.35 -p 50 -z 0 -v 0 -b 0 -qid 10 -cov 10 -E 1 {params.ECOD_N_ITERATIONS} > {output}.log 2>&1'


# annotation
rule annotation:
    input: 
        trigger_search,
        confident_orfs=Path(ORF_PREDICTION_DIR, 'confident_orfs.csv'),
        PCs2proteins=Path(ANNOTATION_DIR, 'PCs2proteins.tsv')
    output: 
        search=Path(OUTPUT_DIR, 'search.tsv'),
        report=Path(OUTPUT_DIR, 'report.tsv'),
        annotation=Path(OUTPUT_DIR, 'annotation.tsv')
    params:
        PHROGS_TABLE=PHROGS_TABLE,
        MGG_PHROGS_TABLE=MGG_PHROGS_TABLE,
        ALAN_TABLE=ALAN_TABLE,
        phrogs_dir=Path(HHSUITE_DIR, 'PHROGS'),
        alan_dir=Path(HHSUITE_DIR, 'ALANDB'),
        pfam_dir=Path(HHSUITE_DIR, 'PFAM'),
        ecod_dir=Path(HHSUITE_DIR, 'ECOD'),
    script: 'dependencies/scripts/annotation.py'


# genbank
rule genbank:
    input:
        annotation=Path(OUTPUT_DIR, 'annotation.tsv'),
        metadata=Path(ORF_PREDICTION_DIR, 'metadata.tsv')
    output: directory(Path(GENBANK_DIR))
    script: 'dependencies/scripts/genbank.py'