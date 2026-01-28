import pandas as pd
import numpy as np

def get_phrog(row):
    """ get phrog cluster """
    if 'phrog_' in str(row['target']):
        return int(row['target'].split('_')[-1])
    else:
        return 0
    
    
def curate_columns(row):
    
    phrog = str(row['phrog'])
    alan_profile = str(row['alan_profile'])
    
    if phrog != '0': return phrog
    elif alan_profile != '0': return alan_profile
    else: return '0'


def combine_function_and_confidence(row):

    function = row['report_function']
    confidence = row['report_confidence']
    
    if function == '-': return function
    else: return f'{confidence:<5}{function}'


def report_phrogs(df, max_evalue=10**-3, nfunc2report=2, verbose=False):

    """
    report a specific number (nfunc2report) of unique functions that got significat (max_evalue) hit(s).
    However, firsty report function significant (not in non-informative list of functions), secondly non-informative functions, lastly unknown functions.
    
    Algorithm:
    1. Filter eval 10**-3
    2. Remove unknown function [REPORT ONLY WHEN NO OTHER FUNCTION, PRIORITY-0]
    3. Remove non-informative functions (lytic tail protein, tail protein, structural protein, virion structural protein, minor tail protein ...) [REPORT ONLY WHEN NO OTHER FUNCTION; PRIORITY-1]
    4. Group by unique functions. For each function report independently max bitscore and max qcov (hits to this function).
    5. Take two functions with highest bitscores.
    6. Report {confidence} {function} in one genbank field, and in seperate field bitscore and qcov.
    7. Report two best PHROG hits seperataly (in total four PHROGS field: 2x function with confidence and 2x params: bitscore and qcov) [PRIORITY-2]
    
    """
    
    # get PC name
    pcid = df['query'].unique()[0]
    
    ### get filters
    noninformative_functions = ['lytic tail protein', 'tail protein', 'structural protein', 'virion structural protein', 'minor tail protein']

    filt_evalue = 'evalue <= @max_evalue'
    get_unknown = '(annot == "unknown function")'
    get_noninformative_functions = '(annot.isin(@noninformative_functions))'
    remove_unknown = '~' + get_unknown
    remove_noninformative = '~' + get_noninformative_functions
    
    informative_query = ' and '.join([remove_unknown, remove_noninformative])
    noninformative_query = get_noninformative_functions

    ### significat only
    df = df.query(filt_evalue)
    
    ### get informative hits
    informative_df = df.query(informative_query)
    informative_df = get_unique_functions_frame(informative_df, function_column='annot') # best hit [max bit score] & highest qcov

    ### noninformative hits
    noninformative_df = df.query(get_noninformative_functions)
    noninformative_df = get_unique_functions_frame(noninformative_df, function_column='annot') # best hit [max bit score] & highest qcov
        
    ### uknown hits
    unknown_df = df.query(get_unknown)
    unknown_df = get_unique_functions_frame(unknown_df, function_column='annot') # best hit [max bit score] & highest qcov

    ### report function
    top_hits_df = pd.concat([informative_df, noninformative_df, unknown_df]).iloc[:nfunc2report].copy()
    
    # Handle empty/short results by padding with no-hit rows
    if len(top_hits_df) < nfunc2report:
        pad = get_no_hit_row(nfunc2report - len(top_hits_df))
        pad['report_function'] = '-'
        # bits/evalue in pad are 0
        top_hits_df = pd.concat([top_hits_df, pad], ignore_index=True)
    
    # Truncate if too many
    top_hits_df = top_hits_df.iloc[:nfunc2report].copy()

    # Labels
    labels = [f'PHROGS{i}' for i in range(1, len(top_hits_df)+1)]
    top_hits_df['query'] = pcid
    top_hits_df['report_label'] = labels
    top_hits_df['report_function'] = top_hits_df['annot']
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
    
    # DO NOT create 'report_params' string here anymore
    # Ensure bits/evalue are numeric
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits']).fillna(0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue']).fillna(0)
    
    return top_hits_df


def report_alan(df, min_prob=0.95, nfunc2report=2, verbose=True):
    """ report a specific number (nfunc2report) of unique functions that got significat (min_prob) hit(s). """
    
    # get PC name
    pcid = df['query'].unique()[0]
    
    # significant & best hits
    df = df.query('prob >= @min_prob')
    df = get_unique_functions_frame(df, function_column='category')
    top_hits_df = df.sort_values('bits', ascending=False).iloc[:nfunc2report].copy()
    
    ### report function
    # prepare columns2report
    if len(top_hits_df) == 0: 
        top_hits_df = get_no_hit_row(2)    
    elif len(top_hits_df) == 1: 
        holder_row = get_no_hit_row(1)
        top_hits_df = pd.concat([top_hits_df, holder_row])
    else: pass

    # report columns
    labels = [f'ALAN{i}' for i in range(1,len(top_hits_df)+1)]
    
    top_hits_df['query'] = [pcid] * len(top_hits_df)
    top_hits_df['report_label'] = labels
    top_hits_df['report_function'] = top_hits_df['category']
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits']).fillna(0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue']).fillna(0)
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
            
    if verbose: display(top_hits_df)

    return top_hits_df

def get_ecod_x(row):
    """
    Extract ECOD X-level from the 'name' field.

    Example ECOD name (generic):
      A|FXXXX|...|X: X-group-id, H: ..., T: ..., F: ...

    This function should return the X-group identifier as a string.
    """
    name = str(row.get("name", ""))
    try:
        parts = name.split("|")
        # parts[3] is usually the hierarchical annotation text
        levels = parts[3].split(": ")
        # "... X: <X-ID>, H: ..."  -> pick the X element
        # levels[2] should contain something like "X: <X-ID>, H"
        x_part = levels[2]          # e.g. "X: 1234, H"
        x_id = x_part.split(",")[0] # "X: 1234"
        x_id = x_id.split()[-1]     # "1234"
        return x_id
    except Exception:
        return "-"
        
def report_ecod(df, nfunc2report=1, verbose=True):
    """
    Report up to nfunc2report ECOD X-groups for a PC.
    """
    pcid = df["query"].unique()[0]

    if df.empty:
        top_hits_df = get_no_hit_row(nfunc2report)
        top_hits_df["report_function"] = "-"
        top_hits_df["report_params"] = "-"
        # ensure numeric
        top_hits_df['bits'] = 0.0
        top_hits_df['evalue'] = 0.0
    else:
        # Extract X level
        if "ecod_x" not in df.columns:
            df = df.copy()
            df["ecod_x"] = df.apply(get_ecod_x, axis=1)

        # Best hit per X group (max bits)
        best_per_x = (
            df.loc[df.groupby("ecod_x")["bits"].idxmax()]
              .sort_values("bits", ascending=False)
              .copy()
        )

        # Keep at most nfunc2report X groups
        best_per_x = best_per_x.iloc[:nfunc2report].copy()

        # Build function string
        def _ecod_function(row):
            name = row["name"]
            try:
                parts = name.split("|")
                F_INDEX = parts[1].strip()
                ecod_levels = parts[3]
                levels = ecod_levels.split(": ")
                T = levels[4].strip(", F")
                F = levels[5].strip()
                return f"X: {row['ecod_x']} T: {T}, F: {F} [{F_INDEX}]"
            except Exception:
                return name

        best_per_x["report_function"] = best_per_x.apply(_ecod_function, axis=1)
        
        top_hits_df = best_per_x

        # Pad with no-hit rows if fewer than requested
        if len(top_hits_df) < nfunc2report:
            pad = get_no_hit_row(nfunc2report - len(top_hits_df))
            pad["report_function"] = "-"
            pad['bits'] = 0.0
            pad['evalue'] = 0.0
            top_hits_df = pd.concat([top_hits_df, pad], ignore_index=True)

    # Labels: ECOD1, ECOD2, ...
    labels = [f"ECOD{i}" for i in range(1, len(top_hits_df) + 1)]
    top_hits_df["query"] = pcid
    top_hits_df["report_label"] = labels
    top_hits_df["report_confidence"] = get_confidence_column(top_hits_df)
    
    # Ensure numeric columns are strictly numeric
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits']).fillna(0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue']).fillna(0)

    if verbose:
        display(top_hits_df)

    return top_hits_df


def report_pfam(df, nfunc2report=1, verbose=True):
    """Report up to nfunc2report PFAM functions (informative first)."""
    pcid = df['query'].unique()[0]

    # Informative first, then DUF, ordered by bits
    informative_df = df.query('~name.str.contains("DUF")').sort_values('bits', ascending=False)
    noninformative_df = df.query('name.str.contains("DUF")').sort_values('bits', ascending=False)
    df_sorted = pd.concat([informative_df, noninformative_df])

    if df_sorted.empty:
        top_hits_df = get_no_hit_row(nfunc2report)
        top_hits_df['report_function'] = '-'
        top_hits_df['bits'] = 0.0
        top_hits_df['evalue'] = 0.0
    else:
        def _pfam_function(row):
            name = row['name']
            try:
                pfamID, func_short, func_detailed = [p.strip() for p in name.split(';')[:3]]
                return f'{func_detailed} [{func_short}] [{pfamID}]'
            except Exception:
                return name

        df_with = df_sorted.copy()
        df_with["report_function"] = df_with.apply(_pfam_function, axis=1)
        
        top_hits_df = df_with.iloc[:nfunc2report].copy()

        if len(top_hits_df) < nfunc2report:
            pad = get_no_hit_row(nfunc2report - len(top_hits_df))
            pad['report_function'] = '-'
            pad['bits'] = 0.0
            pad['evalue'] = 0.0
            top_hits_df = pd.concat([top_hits_df, pad], ignore_index=True)


    labels = [f'PFAM{i}' for i in range(1, len(top_hits_df) + 1)]
    top_hits_df['query'] = pcid
    top_hits_df['report_label'] = labels
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)

    # Ensure numeric columns are strictly numeric
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits']).fillna(0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue']).fillna(0)

    if verbose:
        display(top_hits_df)

    return top_hits_df

def get_no_hit_frames(pcid, nfunc2report):
    """
    Return 4 dataframes (PHROGS, ALAN, PFAM, ECOD) with 'nfunc2report' rows each.
    Labels are numbered (e.g., PHROGS1, PHROGS2) to match hit frames.
    """
    dfs = []
    # Order must match the unpacking in the main script
    db_names = ['PHROGS', 'ALAN', 'PFAM', 'ECOD'] 
    
    for db in db_names:
        # Create n rows
        df = get_no_hit_row(nfunc2report)
        df['query'] = pcid
        df['report_confidence'] = get_confidence_column(df)
        
        # Assign numbered labels: DB1, DB2, ... DBn
        df['report_label'] = [f"{db}{i+1}" for i in range(nfunc2report)]
        
        # Explicitly set numeric columns to float/int 0
        df['bits'] = 0.0
        df['evalue'] = 0.0
        
        dfs.append(df)
        
    return tuple(dfs)

def get_confidence_column(df, eval_intervals=(10**-10, 10**-5, 10**-3, 1), col='evalue', verbose=False):
    # df = pd.DataFrame({'evalue': [10**-10, 10**-7, 10**-3, 10, 10**2, 0]})

    ### conditions
    conditions, choices = [], []
    
    # zero as seperate category (no hits)
    conditions.append(df[col] == 0)
    if verbose: print(f' == 0')
    
    for i, evalue in enumerate(eval_intervals):
        # lower than first number
        if i == 0:
            conditions.append((df[col] <= evalue))
            if verbose: print(f'<= {evalue:.2E}')

        # between numbers
        else:
            lower_evalue = eval_intervals[i-1]
            conditions.append((df[col] > lower_evalue) & (df[col] <= evalue))
            if verbose: print(f'> {lower_evalue:.2E} and <= {evalue:.2E}')
    
    # higher than last number
    conditions.append(df[col] >= eval_intervals[-1])
    if verbose: print(f'> {eval_intervals[-1]}')

    ### choices
    choices = ['*' * i for i in range(len(eval_intervals)).__reversed__()]
    choices = ['-'] + choices[:-1] + ['!', '!']

    # print
    if verbose: print(choices)
    if verbose: np.select(conditions, choices, default='?')

    return np.select(conditions, choices, default='?')


def get_no_hit_row(n, columns_mapper = {'query': 'string', 'target': 'string', 'prob': 'float', \
                                        'pvalue': 'float', 'ident': 'float', 'qcov': 'float', \
                                        'tcov': 'float', 'bits': 'float', 'qstart': 'int', \
                                        'qend': 'int', 'qlength': 'int', 'tstart': 'int', \
                                        'tend': 'int', 'tlength': 'int', 'evalue': 'float', \
                                        'db': 'string', 'name': 'string', 'color': 'string', \
                                        'annot': 'string', 'category': 'string', 'phrog/alan_profile': 'string',
                                        'report_label': 'string', 'report_function': 'string', 'report_params': 'string'}):
    
    """ Give dict of column names and variable types to create 'no hit' row as data frame object """

    values, indicies = [], columns_mapper.keys()
    for key, variable_type in columns_mapper.items():
        if variable_type == 'string': values.append('-')
        else: values.append(0)

    no_hit_row = pd.Series(values, index=indicies).to_frame().T
    no_hit_row = pd.concat([no_hit_row]*n)
    return no_hit_row


def get_unique_functions_frame(df, function_column='annot'):
    """ for each function in data frame:
    - get best hit [highest bitscore]
    - for get highest qcov for a given function from all of the hits
    - return frame of best hits for each function with highest qcov """
    
    # select best hits for each unique function (highest bitscore)
    best_hits_df = df.loc[df.groupby(function_column)['bits'].idxmax()] \
                                                             .sort_values('bits', ascending=False) \
                                                             .copy()
    
    # get highest qcov for each unique function
    best_qcov_df = df.loc[df.groupby(function_column)['qcov'] \
                            .idxmax()][[function_column,'qcov']] \
                            .copy()
    
    # best hits for each unique function & highest qcov for given function
    final_df = best_hits_df.merge(best_qcov_df, on=function_column, how='left', suffixes=('_oryginal', '_best')) \
                           .drop('qcov_oryginal', axis=1) \
                           .rename(columns={'qcov_best': 'qcov'}) \
                           .sort_values('bits', ascending=False)
    
    return final_df