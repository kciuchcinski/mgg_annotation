"""
Generate genbank files of prophages.
"""

# modules
import pandas as pd
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.SeqRecord import SeqRecord


# paths
annotation = snakemake.input.annotation
metadata = snakemake.input.metadata
genbank_dir = snakemake.output[0]

# create
Path(genbank_dir).mkdir(exist_ok=True, parents=True)

# load
annotation_df = pd.read_csv(annotation, sep='\t')
metadata_df = pd.read_csv(metadata, sep='\t')

# curate
annotation_df['contigID'] = annotation_df['proteinID'].str.split(r'_PROTEIN_*', expand=True)[0]


### create genbank
for phageID, group in annotation_df.groupby('contigID'):
    
    # clear features
    features = []

    # phage nucleotide sequence
    seq = metadata_df[metadata_df['contigID'] == phageID]['seq'].values[0]

    # output
    genbank = Path(genbank_dir, f'{phageID}.gb')

    # phage features
    for i, row in group.iterrows():

        start, stop, strand = row['start'], row['stop'], row['strand']
        pc, proteinID, protein = row['PC'], row['proteinID'], row['protein']

        PHROG1_function, PHROG2_function, PHROG1_params, PHROG2_params = row['PHROGS1_function'], row['PHROGS2_function'], row['PHROGS1_params'], row['PHROGS2_params']
        ALAN1_function, ALAN2_function, ALAN1_params, ALAN2_params = row['ALAN1_function'], row['ALAN2_function'], row['ALAN1_params'], row['ALAN2_params']
        PFAM_function, PFAM_params = row['PFAM_function'], row['PFAM_params']
        ECOD_function, ECOD_params = row['ECOD_function'], row['ECOD_params']

        # correct for strand
        if strand == '-': start, stop = stop, start

        f = SeqFeature(FeatureLocation(start-1, stop, strand=int(f'{strand}1')), type='CDS')
        qualifiers = {
                'PHROG1': [PHROG1_function],
                'PHROG2': [PHROG2_function],
                'ALAN1': [ALAN1_function],
                'ALAN2': [ALAN2_function],
                'PFAM': [PFAM_function],
                'ECOD': [ECOD_function],
                'PHROG1_params': [PHROG1_params],
                'PHROG2_params': [PHROG2_params],
                'ALAN1_params': [ALAN1_params],
                'ALAN2_params': [ALAN2_params],
                'PFAM_params': [PFAM_params],
                'ECOD_params': [ECOD_params],
                'proteinID': [proteinID],
                'PC': [pc],
                'seq': [protein],
                'translation': [protein],
                }

        f.qualifiers = qualifiers
        features.append(f)

# This created the duplicated CDS at the end of each genbank file
#     start, stop, strand = row['start'], row['stop'], row['strand']
#     pc, proteinID, protein = row['PC'], row['proteinID'], row['protein']

#     PHROG1_function, PHROG2_function, PHROG1_params, PHROG2_params = row['PHROGS1_function'], row['PHROGS2_function'], row['PHROGS1_params'], row['PHROGS2_params']
#     ALAN1_function, ALAN2_function, ALAN1_params, ALAN2_params = row['ALAN1_function'], row['ALAN2_function'], row['ALAN1_params'], row['ALAN2_params']
#     PFAM_function, PFAM_params = row['PFAM_function'], row['PFAM_params']
#     ECOD_function, ECOD_params = row['ECOD_function'], row['ECOD_params']

#     # correct for strand
#     if strand == '-': start, stop = stop, start

#     f = SeqFeature(FeatureLocation(start-1, stop, strand=int(f'{strand}1')), type='CDS')
#     qualifiers = {
#             'PHROG1': [PHROG1_function],
#             'PHROG2': [PHROG2_function],
#             'ALAN1': [ALAN1_function],
#             'ALAN2': [ALAN2_function],
#             'PFAM': [PFAM_function],
#             'ECOD': [ECOD_function],
#             'PHROG1_params': [PHROG1_params],
#             'PHROG2_params': [PHROG2_params],
#             'ALAN1_params': [ALAN1_params],
#             'ALAN2_params': [ALAN2_params],
#             'PFAM_params': [PFAM_params],
#             'ECOD_params': [ECOD_params],
#             'proteinID': [proteinID],
#             'PC': [pc],
#             'seq': [protein],
#             'translation': [protein],
#             }

#     f.qualifiers = qualifiers
#     features.append(f)

    ### create record
    record = SeqRecord(Seq(seq), id=phageID)
    record.annotations['molecule_type'] = 'DNA'
    record.features = features
    record.name = ''
    record.description = ''
    record.organism = ''

    SeqIO.write(record, genbank, 'genbank')
