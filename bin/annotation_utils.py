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
    pcid = df['query'].unique()[0]
    noninformative_functions = [
        'lytic tail protein', 'tail protein', 'structural protein',
        'virion structural protein', 'minor tail protein'
    ]

    df2 = df.copy()
    df2['evalue'] = pd.to_numeric(df2['evalue'], errors='coerce').fillna(np.inf)
    df2 = df2[df2['evalue'] <= max_evalue]

    if df2.empty:
        top_hits_df = get_no_hit_row(nfunc2report)
        top_hits_df['query'] = pcid
        top_hits_df['report_label'] = [f'PHROGS{i}' for i in range(1, nfunc2report+1)]
        top_hits_df['report_function'] = '-'
        top_hits_df['bits'] = 0.0
        top_hits_df['evalue'] = 0.0
        top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
        return top_hits_df

    unknown_mask = df2['annot'].eq("unknown function")
    noninfo_mask = df2['annot'].isin(noninformative_functions)

    informative_df = get_unique_functions_frame(df2[~unknown_mask & ~noninfo_mask], function_column='annot')
    noninformative_df = get_unique_functions_frame(df2[noninfo_mask], function_column='annot')
    unknown_df = get_unique_functions_frame(df2[unknown_mask], function_column='annot')

    top_hits_df = pd.concat([informative_df, noninformative_df, unknown_df], ignore_index=True).iloc[:nfunc2report].copy()

    if len(top_hits_df) < nfunc2report:
        pad = get_no_hit_row(nfunc2report - len(top_hits_df))
        pad['report_function'] = '-'
        pad['bits'] = 0.0
        pad['evalue'] = 0.0
        top_hits_df = pd.concat([top_hits_df, pad], ignore_index=True)

    top_hits_df = top_hits_df.iloc[:nfunc2report].copy()
    top_hits_df['query'] = pcid
    top_hits_df['report_label'] = [f'PHROGS{i}' for i in range(1, len(top_hits_df)+1)]
    top_hits_df['report_function'] = top_hits_df['annot']
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits'], errors='coerce').fillna(0.0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue'], errors='coerce').fillna(0.0)
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
    return top_hits_df



def report_alan(df, min_prob=0.95, nfunc2report=2, verbose=False):
    pcid = df['query'].unique()[0]
    df2 = df.copy()
    df2['prob'] = pd.to_numeric(df2['prob'], errors='coerce').fillna(0.0)
    df2 = df2[df2['prob'] >= min_prob]

    df2 = get_unique_functions_frame(df2, function_column='category')
    top_hits_df = df2.sort_values('bits', ascending=False).iloc[:nfunc2report].copy()

    if len(top_hits_df) < nfunc2report:
        pad = get_no_hit_row(nfunc2report - len(top_hits_df))
        pad['report_function'] = '-'
        pad['bits'] = 0.0
        pad['evalue'] = 0.0
        top_hits_df = pd.concat([top_hits_df, pad], ignore_index=True)

    top_hits_df['query'] = pcid
    top_hits_df['report_label'] = [f'ALAN{i}' for i in range(1, len(top_hits_df)+1)]
    top_hits_df['report_function'] = top_hits_df.get('category', '-')
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits'], errors='coerce').fillna(0.0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue'], errors='coerce').fillna(0.0)
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
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


def report_pfam(df, nfunc2report=1, verbose=False):
    pcid = df['query'].unique()[0]
    df2 = df.copy()

    # Avoid query('name.str.contains') overhead
    name = df2['name'].astype('string')
    is_duf = name.str.contains("DUF", na=False)

    informative_df = df2[~is_duf].sort_values('bits', ascending=False)
    noninformative_df = df2[is_duf].sort_values('bits', ascending=False)
    df_sorted = pd.concat([informative_df, noninformative_df], ignore_index=True)

    if df_sorted.empty:
        top_hits_df = get_no_hit_row(nfunc2report)
        top_hits_df['query'] = pcid
        top_hits_df['report_label'] = [f'PFAM{i}' for i in range(1, nfunc2report+1)]
        top_hits_df['report_function'] = '-'
        top_hits_df['bits'] = 0.0
        top_hits_df['evalue'] = 0.0
        top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
        return top_hits_df

    def _pfam_function(row):
        name = row['name']
        try:
            pfamID, func_short, func_detailed = [p.strip() for p in str(name).split(';')[:3]]
            return f'{func_detailed} [{func_short}] [{pfamID}]'
        except Exception:
            return str(name)

    top_hits_df = df_sorted.iloc[:nfunc2report].copy()
    if len(top_hits_df) < nfunc2report:
        pad = get_no_hit_row(nfunc2report - len(top_hits_df))
        pad['report_function'] = '-'
        pad['bits'] = 0.0
        pad['evalue'] = 0.0
        top_hits_df = pd.concat([top_hits_df, pad], ignore_index=True)

    top_hits_df['report_function'] = top_hits_df.apply(_pfam_function, axis=1)
    top_hits_df['query'] = pcid
    top_hits_df['report_label'] = [f'PFAM{i}' for i in range(1, len(top_hits_df)+1)]
    top_hits_df['bits'] = pd.to_numeric(top_hits_df['bits'], errors='coerce').fillna(0.0)
    top_hits_df['evalue'] = pd.to_numeric(top_hits_df['evalue'], errors='coerce').fillna(0.0)
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
    return top_hits_df

def get_no_hit_frames(pcid, nfunc2report):
    dfs = []
    db_names = ['PHROGS', 'ALAN', 'PFAM', 'ECOD']
    for db in db_names:
        df = get_no_hit_row(nfunc2report)
        df['query'] = pcid
        df['report_label'] = [f"{db}{i+1}" for i in range(nfunc2report)]
        df['report_function'] = '-'
        df['report_confidence'] = '-'   # no-hit always '-'
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


def get_no_hit_row(n, columns_mapper=None):
    if columns_mapper is None:
        columns_mapper = {
            'query': 'string', 'target': 'string', 'prob': 'float',
            'pvalue': 'float', 'ident': 'float', 'qcov': 'float',
            'tcov': 'float', 'bits': 'float', 'qstart': 'int',
            'qend': 'int', 'qlength': 'int', 'tstart': 'int',
            'tend': 'int', 'tlength': 'int', 'evalue': 'float',
            'db': 'string', 'name': 'string', 'color': 'string',
            'annot': 'string', 'category': 'string', 'phrog/alan_profile': 'string',
            'report_label': 'string', 'report_function': 'string', 'report_params': 'string'
        }

    cols = list(columns_mapper.keys())
    row = [('-' if columns_mapper[c] == 'string' else 0) for c in cols]
    # Much faster than concat([series]*n)
    return pd.DataFrame([row] * n, columns=cols)


def get_unique_functions_frame(df, function_column='annot'):
    """
    For each function:
      - pick best hit by max bits
      - set qcov to max qcov observed for that function
    Avoids extra groupby+merge.
    """
    if df.empty:
        return df.copy()

    df2 = df.copy()
    df2['bits'] = pd.to_numeric(df2['bits'], errors='coerce').fillna(0.0)
    df2['qcov'] = pd.to_numeric(df2['qcov'], errors='coerce').fillna(0.0)

    idx = df2.groupby(function_column, sort=False)['bits'].idxmax()
    best = df2.loc[idx].copy()

    qcov_max = df2.groupby(function_column, sort=False)['qcov'].max()
    best['qcov'] = best[function_column].map(qcov_max)

    return best.sort_values('bits', ascending=False)
