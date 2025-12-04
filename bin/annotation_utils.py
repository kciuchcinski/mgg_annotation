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
    
    # prepare columns2report
    if len(top_hits_df) == 0:
        top_hits_df = get_no_hit_row(2)    
        top_hits_df['report_params'] = ['-', '-']
    elif len(top_hits_df) == 1: 
        holder_row = get_no_hit_row(1)
        holder_row['report_params'] = '-'

        bits, qcov = top_hits_df["bits"].iloc[0], top_hits_df["qcov"].iloc[0]
        top_hits_df['report_params'] = f'bits: {int(bits):>5}  qcov: {qcov:.2f}'
        top_hits_df = pd.concat([top_hits_df, holder_row])
    else: 
        top_hits_df['report_params'] = top_hits_df.apply(lambda row: f'bits: {int(row["bits"]):>5}  qcov: {row["qcov"]:>6.2f}', axis=1)

    # report columns
    labels = [f'PHROGS{i}' for i in range(1,len(top_hits_df)+1)]
    
    top_hits_df['query'] = [pcid] * len(top_hits_df)
    top_hits_df['report_label'] = labels
    top_hits_df['report_function'] = top_hits_df['annot']
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
            
    if verbose: display(top_hits_df)
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
    top_hits_df['report_params'] = top_hits_df.apply(lambda row: f'bits: {int(row["bits"]):>5}  eval: {row["evalue"]:>6.1E}', axis=1)
    top_hits_df['report_confidence'] = get_confidence_column(top_hits_df)
            
    if verbose: display(top_hits_df)

    return top_hits_df



def report_ecod(df, verbose=True):
    """ report best hit [max bitscore] """
    
    # get PC name
    pcid = df['query'].unique()[0]
    
    # sort significant hits
    df = df.sort_values('bits', ascending=False)
    
    # best hit
    top_hit_df = df.iloc[0].copy()
    
    # report function    
    if len(top_hit_df) != 0: 
        name = top_hit_df['name']
        F_INDEX, ecod_levels = name.split('|')[1].strip(), name.split('|')[3]
        T, F = ecod_levels.split(': ')[4].strip(', F'), ecod_levels.split(': ')[5].strip()
        report_function = f'T: {T}, F: {F} [{F_INDEX}]'
        report_params = f'bits: {int(top_hit_df["bits"]):>5}  eval: {top_hit_df["evalue"]:>6.1E}'
    else: # no hit
        top_hits_df = get_no_hit_row(1).iloc[0]
        report_function, report_params = '-', '-'

    # report columns
    top_hit_df['query'] = pcid
    top_hit_df['report_label'] = 'ECOD'
    top_hit_df['report_function'] = report_function
    top_hit_df['report_params'] = report_params
    top_hit_df['report_confidence'] = get_confidence_column(top_hit_df)
            
    if verbose: display(top_hit_df.to_frame().T)
    return top_hit_df.to_frame().T


def report_pfam(df, verbose=True):
    
    # get PC name
    pcid = df['query'].unique()[0]
    
    # sort hits
    informative_df = df.query('~name.str.contains("DUF")').sort_values('bits', ascending=False)
    noninformative_df = df.query('name.str.contains("DUF")').sort_values('bits', ascending=False)
    df = pd.concat([informative_df, noninformative_df])
    
    # best hit & function2report
    if len(df) != 0: 
        top_hit_df = df.iloc[0].copy()
        name = top_hit_df['name']
        pfamID, func_short, func_detailed = name.split(';')[0].strip(), name.split(';')[1].strip(), name.split(';')[2].strip()
        report_function = f'{func_detailed} [{func_short}] [{pfamID}]'
        report_params = f'bits: {int(top_hit_df["bits"]):>5}  eval: {top_hit_df["evalue"]:>6.1E}'
    else: # no hit 
        top_hit_df = get_no_hit_row(1).iloc[0]
        report_function, report_params = '-', '-'
    
    
    # report columns
    top_hit_df['query'] = pcid
    top_hit_df['report_label'] = 'PFAM'
    top_hit_df['report_function'] = report_function
    top_hit_df['report_params'] = report_params
    top_hit_df['report_confidence'] = get_confidence_column(top_hit_df)
            
    if verbose: display(top_hit_df.to_frame().T)
    return top_hit_df.to_frame().T



def get_no_hit_frames(pcid):

    # by defalut: no hit
    phrogs_df, alan_df, pfam_df, ecod_df = get_no_hit_row(2), get_no_hit_row(2), get_no_hit_row(1), get_no_hit_row(1)

    dfs = [phrogs_df, alan_df, pfam_df, ecod_df]
    report_labels = [['PHROGS1', 'PHROGS2'], ['ALAN1', 'ALAN2'], ['PFAM'], ['ECOD']]

    for df, labels in zip(dfs, report_labels):
        df['query'] = pcid
        df['report_confidence'] = get_confidence_column(df)
        df['report_label'] = labels
        
    return phrogs_df, alan_df, pfam_df, ecod_df

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