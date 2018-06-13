import pandas as pd
import os
import os.path
from Bio import SeqIO
from Bio import Entrez

Entrez.email = 'hdutcher@pdx.edu'

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return end1 >= start2 and end2 >= start1

cmsearch = ['target_name', 'target_accession', 'query_name', 'query_accession', \
            'mdl', 'mdl_from', 'mdl_to', 'seq_from', 'seq_to', 'strand', 'trunc', \
            'pass', 'gc', 'bias', 'score', 'E-value', 'inc', 'description of target']

for fn in os.listdir('.'):
    
    if fn.endswith('.txt'):

        print(fn)

        df = pd.read_csv(fn, delim_whitespace=True, header=None, names=cmsearch, \
            skiprows=2, skipfooter=10, engine='python', index_col=False)
    
        # Add columns to reflect minimum and maximum as start/end
        df['start'] = df[['seq_from','seq_to']].min(axis=1)
        df['end'] = df[['seq_from','seq_to']].max(axis=1)

        repeats_reached = False
        seen = []
        hoi = False

        for i, row in df.iterrows():
            
            if not repeats_reached:
            
                if row['target_name'] not in seen:
                    if row['E-value'] <= 0.00001:
                        seen.append(row['target_name'])
                    
                    else:
                        print('Slice results found. Parsing genbank corresponding with \n %s' % row)
                        hoi = True

                        if os.path.isfile(row['target_name'] + '.gb'):
                            foi = open(row['target_name'] + '.gb')
                        else:
                            foi = Entrez.efetch(db='nucleotide', id=row['target_name'], rettype='gbwithparts', retmode='text')
                        
                        with foi as handle:
                            print('parsing genbank...')
                            
                            gb = SeqIO.read(handle, 'genbank')
                            print(gb.description)
                            ol = False
                            for feat in gb.features:
                                if feat.type == 'source':
                                    continue
                                if feat.type == 'gene':
                                    continue
                                else:
                                    try:
                                        if overlap(int(row['start']), int(row['end']), int(feat.location.start), int(feat.location.end)):
                                            print('Overlaps with \n %s' % feat)
                                            ol = True
                                            break

                                    except:
                                        print('Feature skipped due to a problem: %s' % feat)
                                        continue
                            if not ol:
                                print('Hit is intergenic')

                else:
                    repeats_reached = True

            elif not hoi:
                # Find first non-repeating hit
                
                if row['target_name'] not in seen:
                    print('Investigating first non-repeat hit... \n %s' % row)
                    hoi = True
                
                    if os.path.isfile(row['target_name'] + '.gb'):
                        foi = open(row['target_name'] + '.gb')
                    else:
                        foi = Entrez.efetch(db='nucleotide', id=row['target_name'], rettype='gbwithparts', retmode='text')

                    with foi as handle:
                        print('parsing genbank...')
                        
                        gb = SeqIO.read(handle, 'genbank')
                        print(gb.description)
                        ol = False
                        for feat in gb.features:
                            if feat.type == 'gene':
                                continue
                            if feat.type == 'source':
                                continue
                            else:
                                try:
                                    if overlap(int(row['start']), int(row['end']), int(feat.location.start), int(feat.location.end)):
                                        print('Overlaps with \n %s' % feat)
                                        ol = True
                                        break

                                except:
                                    print('Feature skipped due to a problem: %s' % feat)
                                    continue
                        if not ol:
                            print('Hit is intergenic')
                    break

        if not hoi:
            print('No hits of interest found')


