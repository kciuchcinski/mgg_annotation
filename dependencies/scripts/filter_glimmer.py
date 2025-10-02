"""

The code checks if glimmer has ambigous predictions
i.e predictions covering the almost all cds.

The code removes these ambigous cases.

"""

# import modules
import pandas as pd

# paths
glimmer = snakemake.input[0]
filtered = snakemake.output[0]
removed = snakemake.output[1]


# load table
df = pd.read_csv(glimmer)

# reformating columns
df['start'] = pd.to_numeric(df['start'], downcast='integer')
df['stop'] = pd.to_numeric(df['stop'], downcast='integer')
df['diff'] = abs(df['stop'] - df['start'])

# calculating the total cds lenght for every contig
len_cds = {}
for contigID, contig_df in df.groupby(['contigID']):
    total_len = contig_df['stop'].max() - contig_df['start'].min()
    max_cds_len = total_len * 0.9
    len_cds[contigID] = int(max_cds_len)

# filter CDSs that are too long
filt_cds_by_contig = []
for contigID, max_cds_len in len_cds.items():
    filt = (df['contigID'] == contigID)
    filt_cds_by_contig.append(df.loc[filt, 'diff'] >= max_cds_len)

# This part was throwing some errors with Apptainer so I've rewriten it
# In some cases nothing was filtered, so the index len was off

# Initialize a boolean Series with same index as df, all False initially.
filt_too_long = pd.Series(False, index=df.index)

# Loop over each contig, update the Series for CDSs exceeding max length
for contigID, max_cds_len in len_cds.items():
    # Create a boolean Series for this contig
    mask = (df['contigID'] == contigID) & (df['diff'] >= max_cds_len)
    # Update the main Series: set True where condition matches
    filt_too_long.loc[mask] = True


cols = ['start', 'stop', 'strand', 'contigID']

# save filtered data frame
df.loc[~filt_too_long, cols].to_csv(filtered, sep=',', index=False)

# save removed CDSs
df.loc[filt_too_long, cols].to_csv(removed, sep=',', index=False)
