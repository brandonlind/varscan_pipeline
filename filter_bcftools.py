"""Filter output from bcftools."""
from pythonimports import *  # https://github.com/brandonlind/pythonimports/blob/master/pythonimports.py

def filter_snps(df):
    """Mask genotypes if DP<5 or GQ<20, remove loci with >40% missing data.
    
    Parameters
    ----------
    df : pandas.DataFrame
        - rows = loci; output from gatk VariantsToTable
    """
    import numpy as np
    import pandas as pd
    
    # columns for seedlings (megagametophyte) samples
    mg_gtcols = [col for col in df.columns if col.endswith('.GT') and 'mg' in col]
    mg_gqcols = [col for col in df.columns if col.endswith('.GQ') and 'mg' in col]
    mg_dpcols = [col for col in df.columns if col.endswith('.DP') and 'mg' in col]
    
    # columns for parent(s) samples
    p_gtcols = [col for col in df.columns if col.endswith('.GT') and 'mg' not in col]
    p_gqcols = [col for col in df.columns if col.endswith('.GQ') and 'mg' not in col]
    p_dpcols = [col for col in df.columns if col.endswith('.DP') and 'mg' not in col]
    
    # all genotype columns
    all_gtcols = mg_gtcols + p_gtcols
    
    df['frac_missing'] = np.nan
    
    # remove loci that are on multiple lines
    df.index = df['CHROM'] + "-" + df['POS'].astype(str)
    df = df.loc[df.index.value_counts() == 1]
    
    for locus in df.index:
        # for the set of haploid seedlings OR diploid parent(s)
        for gtcols, gqcols, dpcols in [(mg_gtcols, mg_gqcols, mg_dpcols), (p_gtcols, p_gqcols, p_dpcols)]:
#             print(len(gtcols), len(gqcols), len(dpcols))

            # get boolean lists for genotypes passing GQ and DP thresholds
            gqs = df.loc[locus, gqcols] < 20
            gqs.index = gqs.index.str.replace('GQ', 'GT')

            dps = df.loc[locus, dpcols] < 5
            dps.index = dps.index.str.replace('DP', 'GT')

            gts = df.loc[locus, gtcols].isin(['.', './.'])

            # if a locus does not pass thresholds, mask it's genotype
            print('getting good cols')
#             try:
            goodcols = pd.Series(gtcols, index=gtcols)[(gqs) | (dps) | (gts)].tolist()
#             except:
#                 return gtcols, gqs, dps, gts
            if len(goodcols) > 0:
                df.loc[locus, goodcols] = np.nan

        df.loc[locus, 'frac_missing'] = df.loc[locus, all_gtcols].isnull().sum() / len(all_gtcols)
        
    return df


def main(snptable, outdir):
    df = parallel_read(
        snptable,
        lview=lview,
        dview=dview,
        assert_rowcount=False,
        functions=create_fundict(filter_snps)
    )
    
    df = df[(df.index.value_counts() == 1) & (df.frac_missing <=.40)]
    
    outfile = f'{outdir}/%s' % op.basename(snptable).replace('.txt', '_init-filt.txt')
    
    df.to_csv(outfile, index=True, header=True, sep='\t')
    
    print(f'saved filtered snps to\n\t{outfile}')
    
    pass


if __name__ == '__main__':
    thisfile, snptable, outdir, num_engines = sys.argv
    
    lview, dview, cluster_id = start_engines(n=int(num_engines))
    
    main(snptable, outdir)
