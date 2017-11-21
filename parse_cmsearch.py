import pandas as pd
import os
import pickle

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    
    return end1 >= start2 & end2 >= start1

def compare_rows(group):
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
                if winner == j:
                    winners.append(winner)
                    skip.append(i)
                    break
            if last:
                winners.append(i)
    return group.loc[winners].drop_duplicates()

cmsearch = ['target_name', 'target_accession', 'query_name', 'query_accession', \
        'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', \
        'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description of target']

# Load the model name dictionary
with open('cm_dict.p', 'rb') as infile:
    cm_dict = pickle.load(infile)

# Make a master dataframe that contains all the top hits
rdf = pd.DataFrame()
for file in os.listdir('.'):
    if file.endswith('.txt'):
        df = pd.read_csv(file, delim_whitespace=True, header=None, names=cmsearch, \
            skiprows=2, skipfooter=10, engine='python', index_col=False, usecols=cmsearch)

        # Fix the sRNA names if necessary
        spl = lambda x: x.split('_seqs_')[0]
        if file.startswith('RF'):
            df['srna_name'] = df['query_accession'].map(cm_dict)
        else:
            df['srna_name'] = df.query_name.apply(spl)
        
        print(df)

        # Get the top hit only
        fil = df.groupby('target_name')['E-value'].nsmallest(1)

        # Account for weird indexing caused by the groupby operation
        if fil.index.dtype != 'int64':
            fil.index = fil.index.droplevel(level=0)
        keep = df.loc[fil.index]
        print(keep)
        rdf = rdf.append(keep)
        print(rdf)

print(rdf)

# Reset the index to cope with any duplicate indices
rdf = rdf.reset_index()

# Add columns to reflect minimum and maximum as start/end
rdf['start'] = rdf[['seq_from','seq_to']].min(axis=1)
rdf['end'] = rdf[['seq_from','seq_to']].max(axis=1)

# Add column to show sRNA name


print(rdf)
# Remove results that don't cover at least 60% original query length
with open('srna_len.p', 'rb') as infile:
    srna_len = pickle.load(infile)

rdf['srna_orig_len'] = rdf['srna_name'].map(srna_len)

# Calculate the target length
rdf['targ_len'] = rdf['end'].astype(int) - rdf['start'].astype(int) + 1
rdf['perc_coverage'] = rdf['targ_len']/rdf['srna_orig_len']

rdf = rdf.drop(rdf[rdf.perc_coverage < 0.60].index)

# Identify and remove overlaps
rdf_nr = rdf.groupby(['target_name', 'strand'], as_index=False).apply(compare_rows)

# Drop the extra level of indexing gained via the groupby operation
#rdf_nr.index = rdf_nr.index.droplevel(level=0)

# Make and save a presence dictionary
# Create an empty dictionary to hold the genomes for which each sRNA is present
present_dict = {}

srna_grp = rdf_nr.groupby(rdf_nr['srna_name'])

for name, data in srna_grp:
    print(name)
    present = data['target_name'].tolist()
    present_dict[name] = set(acc for acc in present)

# Write presence dictionary to file
with open('rfam_presence_dict.p', 'wb') as outfile:
    pickle.dump(present_dict, outfile)
