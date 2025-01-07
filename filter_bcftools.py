"""Filter output from bcftools."""
from pythonimports import *  # https://github.com/brandonlind/pythonimports/blob/master/pythonimports.py


def remove_dups(df):
    """Remove loci that are on multiple lines."""
    
    locus_dict = df.index.value_counts()
    keep_loci = [locus for locus, count in locus_dict.items() if count == 1]

    df = df[df.index.isin(keep_loci)]

    return df


def mask_snps(df):
    """Mask genotypes if DP<5 or GQ<20, calculate % missing data.
    
    Parameters
    ----------
    df : pandas.DataFrame
        - rows = loci; output from gatk VariantsToTable
    """
    from pythonimports import pbar
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
    
    # column to fill in later
    df['frac_missing'] = np.nan
    
    df.index = df['CHROM'] + "-" + df['POS'].astype(str)
#    df = remove_dups(df)
    
    for locus in pbar(df.index):
        if type(df.loc[locus]) == pd.DataFrame:  # if a locus is on more than one line, skip
            continue
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
            goodcols = pd.Series(gtcols, index=gtcols)[(gqs) | (dps) | (gts)].tolist()
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
        reset_index=False,
        functions=create_fundict(mask_snps)
    )

    print('nrow(df) = ', nrow(df))
    df = remove_dups(df)
    print('nrow(df) = ', nrow(df))
    df = df[df.frac_missing <= 0.40]
    print('nrow(df) = ', nrow(df))
    
    outfile = f'{outdir}/%s' % op.basename(snptable).replace('.txt', '_init-filt.txt')
    
    df.to_csv(outfile, index=True, header=True, sep='\t')
    
    print(f'saved filtered snps to\n\t{outfile}')
    
    pass


if __name__ == '__main__':
    thisfile, snptable, outdir, num_engines = sys.argv
    
    lview, dview, cluster_id = start_engines(n=int(num_engines))
    dview['remove_dups'] = remove_dups
    
    main(snptable, outdir)
