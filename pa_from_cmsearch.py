import sys
import pandas as pd
import pickle
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def pa_matrix(search_results, assemblies_of_interest, refseqinfo):

    df = pickle.load(open(search_results, 'rb'))
    aoi = pickle.load(open(assemblies_of_interest, 'rb'))
    refdict = pickle.load(open(refseqinfo, 'rb'))

    # Create an empty dictionary to hold the genomes for which each sRNA is present
    present_dict = {}

    # Create an empty dataframe to hold the presence/absence matrix
    pa_matrix = pd.DataFrame()
    
    # Group results by sRNA    
    srna_grp = df.groupby(df['query_accession'])

    for name, data in srna_grp:
        print(name)
        pa_dict = {}
        present = data['target_name'].tolist()
        present_dict[name] = set(acc for acc in present)
        for key, value in refdict.items():
            if key in present:
                pa_dict[value] = True
            else:
                if value in pa_dict:
                    continue
                else:
                    pa_dict[value] = False
        # Add that dictionary to the empty data frame
        pa_matrix[name] = pd.Series(pa_dict)

    # Write the presence dictionary to file
    with open('presence_dict.p', 'wb') as outfile:
        pickle.dump(present_dict, outfile)

    # Convert the Trues/Falses to ones and zeros and write this presence/absence matrix to file
    pa_matrix = pa_matrix.astype(int)
    pa_matrix.to_csv('pa_matrix_all_assemblies.csv')
    
    # Cut dataframe to reflect only genomes of interest
    select = pa_matrix.loc[aoi]
    
    # Prepare to convert pa matrix to fasta
    pa_concat = select.apply(lambda row: ''.join(map(str, row)), axis=1)
    
    # Set of records for fasta
    with open('msa_for_GLOOME.fa', 'w') as msa_out:
        msa_records = []
        for i, val in pa_concat.iteritems():
            seq = Seq(str(val))
            record = SeqRecord(seq, id=i)
            msa_records.append(record)
        SeqIO.write(msa_records, msa_out, 'fasta')

if __name__ == '__main__':
    if len(sys.argv) == 4:
         pa_matrix(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print('Usage: pa_from_cmsearch.py nr_search_results.p assemblies_of_interest.p refseqinfo.p')
         sys.exit(0)

