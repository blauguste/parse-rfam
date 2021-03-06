import pickle
import pandas as pd
import sys
import os
import numpy as np

def get_basename(str_in):
    if '_v' in str_in:
        orig = str_in.split('_v')[0]
        return orig
    else:
        return str_in

def get_acc_num(str_in):
    acc_num = str_in.split('|')[1]
    return acc_num

def get_sp_name(str_in):
    split_sp_name = str(str_in).split(' ')
    sp_name = split_sp_name[0] + ' ' + split_sp_name[1]
    return sp_name

def get_shortname(str_in):
    split = str(str_in).split(' ')
    short = split[0][0] + '_' + split[1]
    return short

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return end1 >= start2 and end2 >= start1

def compare_rows(passed_group):
    global redun
    group = passed_group.sort_values(by='E-value', ascending=True)
    print(group)
    winners = []
    skip = []
    if len(group) == 1:
        return group
    for i in group.index:
        if i in skip:
            continue
        for j in group.index:
            last = j == group.index[-1]
            istart = group.loc[i, 'start']
            iend = group.loc[i, 'end']
            jstart = group.loc[j, 'start']
            jend = group.loc[j, 'end']
            if overlap(istart, iend, jstart, jend):
                winner = group.loc[[i, j], 'E-value'].idxmin()
                loser = group.loc[[i, j], 'E-value'].idxmax()
                redun.append((group.loc[winner].query_accession, group.loc[loser].query_accession))
                if winner == j:
                    winners.append(winner)
                    skip.append(i)
                    break
            if last:
                winners.append(i)
    return group.loc[winners].drop_duplicates()

def nr_cmsearch(outfile):

    global redun
    redun = []

    cmsearch = ['target_name', 'target_accession', 'query_name', 'query_accession', \
            'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', \
            'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description of target']

    with open(outfile, 'ab') as ofn:
        writer = pd.ExcelWriter(ofn)

        rdf = pd.DataFrame()

        for fn in os.listdir('.'):
            if fn.endswith(('.tab', '.txt')):
                df = pd.read_csv(fn, delim_whitespace=True, header=None, names=cmsearch, \
                    skiprows=2, skipfooter=10, engine='python', index_col=False)

                # Get top hit only
                fil = df.groupby('target_name')['E-value'].nsmallest(1)
                
                if isinstance(fil.index, pd.core.index.MultiIndex):
                    fil.index = fil.index.droplevel()
                                
                tophit = df.loc[fil.index]

                rdf = rdf.append(tophit)

        # Reset the index to cope with any duplicate indices
        rdf = rdf.reset_index()

        # Drop column named 'index' (not the real index)
        rdf = rdf.drop(['index'], axis=1)
        
        # Add columns to reflect minimum and maximum as start/end
        rdf['start'] = rdf[['seq_from','seq_to']].min(axis=1)
        rdf['end'] = rdf[['seq_from','seq_to']].max(axis=1)

        # Remove results that don't cover at least 60% original query length
        with open('srna_len.p', 'rb') as infile:
            srna_len = pickle.load(infile)

        # Put identifiers all in same column
        rdf = rdf.replace({'query_accession': {'-': np.nan}})
        rdf.query_accession = rdf.query_accession.fillna(rdf.query_name)
        rdf['basename'] = rdf['query_accession'].apply(get_basename)

        rdf['srna_orig_len'] = rdf['basename'].map(srna_len)
        rdf['targ_len'] = rdf['end'].astype(int) - rdf['start'].astype(int) + 1
        rdf['perc_coverage'] = rdf['targ_len']/rdf['srna_orig_len']

        rdf = rdf.drop(rdf[rdf.perc_coverage < 0.60].index)

        rdf.to_pickle('full_tophit_results.p')

        # Identify and remove overlaps
        rdf_nr = rdf.groupby(['target_name', 'strand'], as_index=False).apply(compare_rows)

        rdf_nr.index = rdf_nr.index.droplevel()

        # Save the nonredundant results for later use
        rdf_nr.to_pickle('nr_cmsearch_results.p')
        
        with open('redundancy_tuples.p', 'wb') as outfile:
            pickle.dump(redun, outfile)

        # Find the difference between the original results and the nr results
        overlaps = rdf[~rdf.index.isin(rdf_nr.index)]

        rdf_nr.to_excel(writer, sheet_name='nr_cmsearch')
        overlaps.to_excel(writer, sheet_name='removed_due_to_overlap')
        writer.save()

if __name__ == '__main__':
    if len(sys.argv) == 2:
         nr_cmsearch(sys.argv[1])
    else:
         print('Usage: nr_cmsearch.py outfile.xlsx')
         sys.exit(0)
