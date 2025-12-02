"""Filter output from bcftools.

TODO
----
- keep multiallelic if REF not in either ALT, and len(ALT.split(',')) == 2
"""
from pythonimports import *  # https://github.com/brandonlind/pythonimports/blob/master/pythonimports.py

def calc_AF(df):
    """Calculate AF and MAF for each locus."""
    gtcols = [col for col in df.columns if col.endswith('.GT')]
    for locus in pbar(df.index):
        alt = df.loc[locus, 'ALT']
        if alt == '.':  # NO_VARIATION, AF set in mask_snps function below
            continue
        samp_gts = df.loc[locus, gtcols]
        gts = samp_gts[samp_gts.notnull()]

        assert gts.str.count(r'\.').sum() == 0  # assumes missing gts are np.nan

        alt_count = gts.str.count(alt).sum()  # eg if alt=A and gt=A/A, count=2

        AF = alt_count / (2 * len(gts))

        if AF == 0:
            df.loc[locus, 'TYPE'] = 'NO_VARIATION'

        df.loc[locus, 'AF'] = AF

    df['MAF'] = df.AF.apply(lambda af: af if af <= 0.5 else 1 - af)
    return df


def remove_dups(df):
    """Remove loci that are on multiple lines.

    Notes
    -----
    - remove_dups is necessary when doing parallel_read because a locus may exist in unique jobs sent to parallel engines.
    - otherwise filtered out if df.loc[locus] is a dataframe instead of a series in mask_snps function below
    """

    locus_dict = df.index.value_counts()
    keep_loci = [locus for locus, count in locus_dict.items() if count == 1]

    df = df[df.index.isin(keep_loci)]

    return df


def mask_snps(df, mindepth=5, min_gq=20):
    """Mask genotypes if DP<min_depth or GQ<min_gq, calculate % missing data.

    Parameters
    ----------
    df : pandas.DataFrame
        - rows = loci; output from gatk VariantsToTable
    """
    from pythonimports import pbar
    import numpy as np
    import pandas as pd

    # columns for samples - assumes sample columns are grouped (ie groups are ordered by sample)
    gtcols = [col for col in df.columns if col.endswith('.GT')]
    gqcols = [col for col in df.columns if col.endswith('.GQ')]
    dpcols = [col for col in df.columns if col.endswith('.DP')]
    adcols = [col for col in df.columns if col.endswith('.AD')]

    # column to fill in later
    df['frac_missing'] = np.nan

    df.index = df['CHROM'] + "-" + df['POS'].astype(str)

    for locus in pbar(df.index):
        if type(df.loc[locus]) == pd.DataFrame:  # if a locus is on more than one line, skip
            continue

        if df.loc[locus, 'TYPE'] in ['SNP', 'INDEL']:
            # get boolean lists for genotypes passing GQ and DP thresholds
            gqs = df.loc[locus, gqcols] < min_gq
            gqs.index = gqs.index.str.replace('GQ', 'GT')

            dps = df.loc[locus, dpcols] < mindepth
            dps.index = dps.index.str.replace('DP', 'GT')

            gts = df.loc[locus, gtcols].isin(['.', './.'])

            # if a locus does not pass thresholds, mask it's genotype
            badcols = pd.Series(gtcols, index=gtcols)[(gqs) | (dps) | (gts)].tolist()
            if len(badcols) > 0:
                df.loc[locus, badcols] = np.nan

            # fix HET counts after masking
            ref, alt = df.loc[locus, ['REF', 'ALT']]
            het1 = f'{ref}/{alt}'
            het2 = f'{alt}/{ref}'
            df.loc[locus, 'HET'] = df.loc[locus, gtcols].str.count(het1).sum() + df.loc[locus, gtcols].str.count(het2).sum()

        elif df.loc[locus, 'TYPE'] == 'NO_VARIATION':  # invariant sites don't output GQs
            # get boolean lists for genotypes passing GQ and DP thresholds

            dps = df.loc[locus, dpcols] < mindepth
            dps.index = dps.index.str.replace('DP', 'GT')

            gts = df.loc[locus, gtcols].isin(['.', './.'])

            # if a locus does not pass thresholds, mask it's genotype
            badcols = pd.Series(gtcols, index=gtcols)[(dps) | (gts)].tolist()
            if len(badcols) > 0:
                df.loc[locus, badcols] = np.nan

            # fix from NAN
            df.loc[locus, 'AF'] = 0
        else:
            raise Exception('unaccounted TYPE: %s' % df.loc[locus, 'TYPE'])

        # fix AD by counting depths for unmasked genotypes
        gts = df.loc[locus, gtcols]
        ads = df.loc[locus, adcols]
        if df.loc[locus, 'TYPE'] == 'NO_VARIATION':
            ref_depth = ads[ads.index[gts.notnull()]].astype(int).sum()
            alt_depth = 0
        else:
            ref_depth = ads[ads.index[gts.notnull()]].apply(lambda x: x.split(',')[0]).astype(int).sum()
            alt_depth = ads[ads.index[gts.notnull()]].apply(lambda x: x.split(',')[1]).astype(int).sum()
        df.loc[locus, 'AD'] = f'{ref_depth},{alt_depth}'

        # calc frac_missing
        df.loc[locus, 'frac_missing'] = df.loc[locus, gtcols].isnull().sum() / len(gtcols)

    return df


def mark_nas(df):
    gtcols = [col for col in df.columns if col.endswith('.GT')]
    for gtcol in gtcols:
        samp = gtcol.removesuffix('.GT')
        sampcols = [col for col in df.columns if col.startswith(f'{samp}.')]
        nas = df[gtcol].isnull()
        df.loc[nas, sampcols] = np.nan
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

    print('starting with: nrow(df) = ', nrow(df))
    df = remove_dups(df)
    print('after remove_dups: nrow(df) = ', nrow(df))
    df = df[df.frac_missing <= 0.40]
    print('after frac_missing: nrow(df) = ', nrow(df))
    df = calc_AF(df)
    print('after calc_AF: nrow(df) = ', nrow(df))
    df = mark_nas(df)
    print('after mark_nas: nrow(df) = ', nrow(df))

    starting_cols = ['CHROM', 'POS', 'REF', 'ALT', 'AF', 'MAF', 'frac_missing', 'QUAL', 'TYPE', 'FILTER', 'HET', 'AD']
    df = df[starting_cols + [col for col in df.columns if col not in starting_cols]]

    outfile = f'{outdir}/%s' % op.basename(snptable).replace('.txt', '_init-filt.txt.gz')

    df.to_csv(outfile, index=True, header=True, sep='\t')

    print(f'saved filtered snps to\n\t{outfile}')

    pass


if __name__ == '__main__':
    thisfile, snptable, outdir, num_engines = sys.argv

    lview, dview, cluster_id = start_engines(n=int(num_engines))

    main(snptable, outdir)
