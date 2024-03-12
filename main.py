"""


    """

import pandas as pd
from bgen_reader import open_bgen
from path import Path
from bgen_reader import read_bgen , allele_expectation , compute_dosage
import numpy as np

class MfiCols :
    maf = 5
    info = 7

mc = MfiCols()

class Dirs :
    mahdi_ukb_dir = Path('/disk/genetics/ukb/mahdimir')
    proj_dir = mahdi_ukb_dir / 'imputed_genotype_corr_Tammys_analysis_replication'
    snps_2_keep = proj_dir / 'snps_2_keep'
    plink_out = proj_dir / 'plink_out_sibs'
    out_tammy_way = proj_dir / 'out_Tammy_s_way'
    out = proj_dir / 'out'
    med_g = Path(
            '/var/genetics/ws/mahdimir/med/imputed-genotype-sibling-task-240311')
    med_gpro = Path(
            '/var/genetics/ws/mahdimir/DropBox/2-GPro/1-med-GPro/imputed-genotype-sibling-task-240311')

dyr = Dirs()

def filter_snps() :
    """
    filter out 1000 snps with maf > 0.01 & ( .3 < info <= .31 ) or snps with maf > 0.01 & ( .99 < info <= 1 )
    """

    ##
    snps_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/unpacked_ukb_imp_mfi_tgz/ukb_mfi_chr22_v3.txt'

    df = pd.read_csv(snps_fp , sep = '\s+' , header = None)

    ##
    # filter snps with maf > 0.01 & ( .3 < info <= .31 )

    # keep snps with maf > 1%
    msk = df[mc.maf].gt(.01)

    # keep snps with ( .3 < info <= .31 )
    msk &= (df[mc.info].gt(.3)) & (df[mc.info].le(.31))

    dfa = df[msk]

    if len(dfa) > 1000 :
        dfa = dfa.sample(1000)

    dfa.to_csv(dyr.med_gpro / 'snps_30.txt' ,
               index = False ,
               sep = '\t' ,
               header = False)

    ##
    # keep snps with maf > 1%
    msk = df[mc.maf].gt(.01)

    # keep snps with .99 < info <= 1
    msk &= (df[mc.info].gt(.99)) & (df[mc.info].le(1))

    dfa = df[msk]

    if len(dfa) > 1000 :
        dfa = dfa.sample(1000)

    dfa.to_csv(dyr.med_gpro + '/snps_99.txt' ,
               index = False ,
               sep = '\t' ,
               header = False)

def filter_sibs() :
    """ """
    ##
    rel_fp = '/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0'

    df = pd.read_csv(rel_fp , sep = '\s+' , dtype = 'string')

    # keep only full sibs
    msk = df['InfType'].eq('FS')
    dfa = df[msk]

    df1 = dfa[['FID1' , 'ID1']]
    df2 = dfa[['FID2' , 'ID2']]

    df1.columns = ['FID' , 'ID']
    df2.columns = ['FID' , 'ID']

    dfb = pd.concat([df1 , df2])
    dfb = dfb.iloc[: , :2]

    dfb.to_csv(dyr.med_gpro + '/sibs.txt' ,
               index = False ,
               header = False ,
               sep = '\t')

def filter_parent_offspring() :
    """ """
    ##
    rel_fp = '/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0'

    df = pd.read_csv(rel_fp , sep = '\s+' , dtype = 'string')

    ##
    # keep only full sibs
    msk = df['InfType'].eq('PO')
    dfa = df[msk]

    df1 = dfa[['FID1' , 'ID1']]
    df2 = dfa[['FID2' , 'ID2']]

    df1.columns = ['FID' , 'ID']
    df2.columns = ['FID' , 'ID']

    dfb = pd.concat([df1 , df2])
    dfb = dfb.iloc[: , :2]

    dfb.to_csv(dyr.med_gpro + '/parent_offspring.txt' ,
               index = False ,
               header = False ,
               sep = '\t')

def compute_corr_dosages() :
    """"""

    ##
    sibs_fp = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0"
    df_sibs = pd.read_csv(sibs_fp , sep = '\s+' , dtype = 'string')

    ##
    msk = df_sibs['InfType'].eq('FS')

    df_fs = df_sibs[msk]

    ##
    raw_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bgen_30.raw'
    df_raw = pd.read_csv(raw_fp , sep = '\s+')

    ##
    df_raw = df_raw.drop(columns = ['FID' , 'PAT' , 'MAT' , 'SEX' ,
                                    'PHENOTYPE'])

    ##
    df_sib1 = pd.merge(df_fs[['ID1']] ,
                       df_raw ,
                       left_on = 'ID1' ,
                       right_on = 'IID' ,
                       how = 'left')
    df_sib2 = pd.merge(df_fs[['ID2']] ,
                       df_raw ,
                       left_on = 'ID2' ,
                       right_on = 'IID' ,
                       how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['IID'])
    df_sib2 = df_sib2.drop(columns = ['IID'])

    ##
    df_sib1 = df_sib1.rename(columns = {
            'ID1' : 'IID'
            })
    df_sib2 = df_sib2.rename(columns = {
            'ID2' : 'IID'
            })

    ##
    df_sib1 = df_sib1.astype('float')
    df_sib2 = df_sib2.astype('float')

    ##
    df_cors = df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##
    df_cors = df_cors.drop(index = 'IID')

    ##
    df_cors.to_csv(dyr.out_tammy_way / 'sib_corr_dosages_30.csv')

    ##

    ##
    sibs_fp = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0"
    df_sibs = pd.read_csv(sibs_fp , sep = '\s+' , dtype = 'string')

    ##
    msk = df_sibs['InfType'].eq('FS')

    df_fs = df_sibs[msk]

    ##
    raw_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bgen_99.raw'
    df_raw = pd.read_csv(raw_fp , sep = '\s+' , dtype = 'string')

    ##
    df_raw = df_raw.drop(columns = ['FID' , 'PAT' , 'MAT' , 'SEX' ,
                                    'PHENOTYPE'])

    ##
    df_sib1 = pd.merge(df_fs[['ID1']] ,
                       df_raw ,
                       left_on = 'ID1' ,
                       right_on = 'IID' ,
                       how = 'left')
    df_sib2 = pd.merge(df_fs[['ID2']] ,
                       df_raw ,
                       left_on = 'ID2' ,
                       right_on = 'IID' ,
                       how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['IID'])
    df_sib2 = df_sib2.drop(columns = ['IID'])

    ##
    df_sib1 = df_sib1.rename(columns = {
            'ID1' : 'IID'
            })
    df_sib2 = df_sib2.rename(columns = {
            'ID2' : 'IID'
            })

    ##
    df_sib1 = df_sib1.astype('float')
    df_sib2 = df_sib2.astype('float')

    ##
    df_cors = df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##
    df_cors = df_cors.drop(index = 'IID')

    ##
    df_cors.to_csv(dyr.out_tammy_way / 'sib_corr_dosages_99.csv')

    ##

def compute_corr_hard_calls() :
    """
    compute correlation between sibs using hard calls
    :return:
    """

    ##
    sibs_fp = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0"
    df_sibs = pd.read_csv(sibs_fp , sep = '\s+' , dtype = 'string')

    ##
    msk = df_sibs['InfType'].eq('FS')

    df_fs = df_sibs[msk]

    ##
    bed_fp = "/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bed/bgen_30.bed"

    bed = Bed(bed_fp , count_A1 = True)

    ##
    snps = bed.read()

    ##
    x = snps.val

    ##
    df_bed_0 = pd.DataFrame(x , columns = bed.sid)

    ##
    fam_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bed/bgen_30.fam'
    df_fam = pd.read_csv(fam_fp ,
                         sep = '\s+' ,
                         dtype = 'string' ,
                         header = None)

    ##
    df_bed = pd.concat([df_fam[[1]] , df_bed_0] , axis = 1)

    ##

    ##
    df_bed = df_bed.rename(columns = {
            1 : 'ID'
            })

    ##
    df_sib1 = pd.merge(df_fs[['ID1']] ,
                       df_bed ,
                       left_on = 'ID1' ,
                       right_on = 'ID' ,
                       how = 'left')

    ##
    df_sib2 = pd.merge(df_fs[['ID2']] ,
                       df_bed ,
                       left_on = 'ID2' ,
                       right_on = 'ID' ,
                       how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['ID' , 'ID1'])
    df_sib2 = df_sib2.drop(columns = ['ID' , "ID2"])

    ##
    df_sib1 = df_sib1.astype('float')
    df_sib2 = df_sib2.astype('float')

    ##
    df_cors = df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##
    df_cors.to_csv(dyr.out_tammy_way / 'sib_corr_hard_calls_30.csv')

    ##

    ##
    sibs_fp = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0"
    df_sibs = pd.read_csv(sibs_fp , sep = '\s+' , dtype = 'string')

    ##
    msk = df_sibs['InfType'].eq('FS')

    df_fs = df_sibs[msk]

    ##
    bed_fp = "/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bed/bgen_99.bed"

    bed = Bed(bed_fp , count_A1 = True)

    ##
    snps = bed.read()

    ##
    x = snps.val

    ##
    df_bed_0 = pd.DataFrame(x , columns = bed.sid)

    ##
    fam_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bed/bgen_99.fam'
    df_fam = pd.read_csv(fam_fp ,
                         sep = '\s+' ,
                         dtype = 'string' ,
                         header = None)

    ##
    df_bed = pd.concat([df_fam[[1]] , df_bed_0] , axis = 1)

    ##
    df_bed = df_bed.rename(columns = {
            1 : 'ID'
            })

    ##
    df_sib1 = pd.merge(df_fs[['ID1']] ,
                       df_bed ,
                       left_on = 'ID1' ,
                       right_on = 'ID' ,
                       how = 'left')

    ##
    df_sib2 = pd.merge(df_fs[['ID2']] ,
                       df_bed ,
                       left_on = 'ID2' ,
                       right_on = 'ID' ,
                       how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['ID' , 'ID1'])
    df_sib2 = df_sib2.drop(columns = ['ID' , "ID2"])

    ##
    df_sib1 = df_sib1.astype('float')
    df_sib2 = df_sib2.astype('float')

    ##
    df_cors = df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##
    df_cors.to_csv(dyr.out_tammy_way / 'sib_corr_hard_calls_99.csv')

    ##

    ##

    ##
    fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'
    from bgen_reader import open_bgen

    bgen = open_bgen(fp , verbose = False)

    ##
    fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_3/bgen_30.bgen'
    bgen = open_bgen(fp , verbose = False)

    ##
    print(bgen.rsids)

    ##
    len(bgen.rsids)

    ##

    ##
    s1 = df_sib1.iloc[: , [0]]
    s2 = df_sib2.iloc[: , [0]]

    # get pairwise correlation
    s1.corrwith(s2 , method = 'pearson')

    ##
    x = df_sib1.iloc[: , [0]].to_numpy()
    y = df_sib2.iloc[: , [0]].to_numpy()

    ##
    from scipy.stats import pearsonr

    pearsonr(df_sib1.iloc[: , 1] , df_sib2.iloc[: , 1])

    ##
    df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##

    ##
    df = pd.DataFrame(np.random.random((5 , 5)) ,
                      columns = ['gene_' + chr(i + ord('a')) for i in range(5)])
    df1 = pd.DataFrame(np.random.random((5 , 5)) ,
                       columns = ['gene_' + chr(i + ord('a')) for i in
                                  range(5)])
    print(df)
    pearsonr(df.loc[: , 'gene_a'] , df1.loc[: , 'gene_a'])

    ##
    import itertools

    correlations = {}
    columns = df.columns.tolist()

    for col_a , col_b in itertools.combinations(columns , 2) :
        correlations[col_a + '__' + col_b] = pearsonr(df.loc[: , col_a] ,
                                                      df.loc[: , col_b])

    ##

    ##
    df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##
    df_sib1[['rs535733662_G']].corrwith(df_sib2[['rs535733662_G']] ,
                                        method = 'pearson')

    ##

    df_cors = df_sib1.corrwith(df_sib2 , method = 'pearson')

    ##
    df_sib1.equals(df_sib2)

    ##

    ##

    ##

    ##
    df_sib1 = df_raw[df_raw['IID'].isin(df_fs['ID1'])]
    df_sib2 = df_raw[df_raw['IID'].isin(df_fs['ID2'])]

    ##

    ##

    ##

##


##

def compare_bgen_with_raw() :
    """ """

    ##
    sibs_fp = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0"
    df_sibs = pd.read_csv(sibs_fp , sep = '\s+' , dtype = 'string')

    ##
    msk = df_sibs['InfType'].eq('FS')

    df_fs = df_sibs[msk]

    ##
    bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'

    # reading bgen files in two different objects with different methods and attrs
    bg_opn = open_bgen(bgen_fp , verbose = False)
    bg_read = read_bgen(bgen_fp , verbose = False)

    # samples
    bg_opn_s = list(bg_opn.samples)

    ##
    # making a data frame to gather dosages from all variants
    df_d = pd.DataFrame({
            'IID' : bg_opn_s
            })

    def ret_dosages(bg_read , variant_idx , alt_allele_idx = 1) :
        e = allele_expectation(bg_read , variant_idx)
        d = compute_dosage(e , alt = alt_allele_idx)
        return d

    for v_idx in range(bg_opn.nvariants) :
        d = ret_dosages(bg_read , v_idx)
        df_d[v_idx] = d

    ##
    df_d['IID'] = df_d['IID'].str.split('_').str[1]

    ##
    bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'

    # reading bgen files in two different objects with different methods and attrs
    bg_opn = open_bgen(bgen_fp , verbose = False)

    ##
    def get_hard_call(bg_read , variant_idx) :
        gts = bg_read['genotype'][variant_idx].compute()
        dfp = pd.DataFrame(gts['probs'])
        sr = dfp.idxmax(axis = 1)
        sr = sr - 2
        sr = sr.abs()
        return sr

    ##
    df_hc = pd.DataFrame({
            'IID' : bg_opn_s
            })

    for v_idx in range(bg_opn.nvariants) :
        df = get_hard_call(bg_read , v_idx)
        df_hc[v_idx] = df

    ##
    df_hc['IID'] = df_hc['IID'].str.split('_').str[1]

    ##
    df_hc.isna().sum().sum()

    ##
    # check whether they are the same except on nan values
    df_hc.equals(df_bed)

    # we have many nan vals

    ##
    bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'

    # reading bgen files in two different objects with different methods and attrs
    bg_opn = open_bgen(bgen_fp , verbose = False)

    ##
    x = bg_opn.read()

    ##
    type(x)

    ##
    x.nbytes / 1024 ** 2

    ##  # get the results firtst then try to remove for the sake of speed

##


def get_sib_pairs_ids() :
    """ """

    ##
    rel_fp = '/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0'
    df = pd.read_csv(rel_fp , sep = '\s+' , dtype = 'string')

    msk = df['InfType'].eq('FS')

    df = df[msk]

    df = df[['ID1' , 'ID2']]

    ##
    return df

##
def ret_dosages(bg_read , variant_idx , alt_allele_idx = 1) :
    e = allele_expectation(bg_read , variant_idx)
    d = compute_dosage(e , alt = alt_allele_idx)
    return d

##
def make_df_of_iids_from_bgen_open_obj(bg_opn) :
    """ """

    ##
    bg_opn_s = list(bg_opn.samples)
    df = pd.DataFrame({
            'IID' : bg_opn_s
            })
    df['IID'] = df['IID'].str.split('_').str[1]

    ##
    return df

##
def get_dosages_of_all_vars_from_bgen_using_for(bgen_fp) :
    """ """
    # this is not an efficient way to do this
    if False :
        pass

        ##
        bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_99.bgen'

    ##
    bg_opn = open_bgen(bgen_fp)
    bg_read = read_bgen(bgen_fp)

    ##
    df2 = make_df_of_iids_from_bgen_open_obj(bg_opn)

    ##
    for v_idx in range(5) :
        d = ret_dosages(bg_read , v_idx)
        df2[v_idx] = d

    ##
    return df

##


##
def get_hard_call_of_a_variant(bg_read , variant_idx) :
    """ """

    if False :
        pass

        ##
        bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_99.bgen'

    ##
    gts = bg_read['genotype'][variant_idx].compute()
    dfp = pd.DataFrame(gts['probs'])
    sr = dfp.idxmax(axis = 1)
    sr = sr - 2
    sr = sr.abs()
    return sr

##
def get_hard_calls_of_all_vars_from_bgen(bgen_fp) :
    """ """

    if False :
        pass

        ##
        bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'

    ##
    bg_opn = open_bgen(bgen_fp)
    bg_read = read_bgen(bgen_fp)

    ##
    df = make_df_of_iids_from_bgen_open_obj(bg_opn)

    ##
    for v_idx in range(bg_opn.nvariants) :
        sr = get_hard_call_of_a_variant(bg_read , v_idx)
        df[v_idx] = sr

    ##
    return df

##
def get_dosages_and_hardcall_data_of_all_vars_from_bgen(bgen_fp) :
    """ """

    if False :
        pass

        ##
        bgen_fp = '/var/genetics/ws/mahdimir/DropBox/2-GPro/1-med-GPro/imputed-genotype-sibling-task-240311/plink_out/snps_99.bgen'

    ##
    bg_opn = open_bgen(bgen_fp)

    ##
    df_id = make_df_of_iids_from_bgen_open_obj(bg_opn)

    ##
    nd_p = bg_opn.read()

    ##
    nd_d = nd_p[: , : , 1] + 2 * nd_p[: , : , 2]

    ##
    df1 = pd.DataFrame(nd_d)
    del nd_d

    ##
    df_d = pd.concat([df_id , df1] , axis = 1)

    ##
    del df1

    ##
    _fp = dyr.out / 'df_d.parquet'
    df_d.to_parquet(_fp)

    ##
    del df_d

    ##
    nd_h = np.argmax(nd_p , axis = 2)

    ##
    df1 = pd.DataFrame(nd_h)
    del nd_h

    ##
    df_h = pd.concat([df_id , df1] , axis = 1)

    ##
    _fp = dyr.out / 'df_h.parquet'
    df_h.to_parquet(_fp)

    ##
    _fpd = dyr.out / 'df_d.parquet'
    _fph = dyr.out / 'df_h.parquet'

    df_d = pd.read_parquet(_fpd)
    df_h = pd.read_parquet(_fph)

    ##
    dfs = [df_d , df_h]

    ##
    methods = {
            0 : 'dosages' ,
            1 : 'hard_calls'
            }

    for df , mtd in zip(dfs , methods.values()) :
        save_corr_of_sib_pairs_with_gts_df(99 , mtd , df)

    ##

    ##
    return df_d , df_h

##
def make_sib_pairs_gts(df_gts , df_sibs) :
    """ """
    if False :
        pass

        ##
        bgen_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'

        ##
        df_gts = get_hard_calls_of_all_vars_from_bgen(bgen_fp)

    ##
    df_sib1 = pd.merge(df_sibs[['ID1']] ,
                       df_gts ,
                       left_on = 'ID1' ,
                       right_on = 'IID' ,
                       how = 'left')
    df_sib2 = pd.merge(df_sibs[['ID2']] ,
                       df_gts ,
                       left_on = 'ID2' ,
                       right_on = 'IID' ,
                       how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['ID1'])
    df_sib2 = df_sib2.drop(columns = ['ID2'])

    ##
    return df_sib1 , df_sib2

##
def save_corr_of_sib_pairs(info_score) :
    """ """
    if False :
        pass

        ##
        info_score = 99

    ##
    bg_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_{}.bgen'
    bg_fp = bg_fp.format(info_score)
    bg_fp

    ##
    _f = get_dosages_and_hardcall_data_of_all_vars_from_bgen
    dfs = _f(bg_fp)

    ##
    methods = {
            0 : 'dosages' ,
            1 : 'hard_calls'
            }

    for df , mtd in zip(dfs , methods.values()) :
        save_corr_of_sib_pairs_with_gts_df(info_score , mtd , df)

    ##

    ##

def save_corr_of_sib_pairs_with_gts_df(info_score , gts_method , df_gts) :
    """ """

    ##
    df_sibs = get_sib_pairs_ids()

    ##
    df_sib1 , df_sib2 = make_sib_pairs_gts(df_gts , df_sibs)

    ##
    gts1 = df_sib1.iloc[: , 1 :]
    gts2 = df_sib2.iloc[: , 1 :]

    ##
    df_cors = gts1.corrwith(gts2 , method = 'pearson')

    ##
    out_fp = dyr.out / f'sib_corr_{gts_method}_{info_score}.csv'
    df_cors.to_csv(out_fp , index = False , header = False)
    print(out_fp)

##
def main() :
    pass

    ##
    info_scores = {
            0 : 30 ,
            1 : 99 ,
            }

    ##
    for info in info_scores.values() :
        pass

        ##
        save_corr_of_sib_pairs(info)

        ##  ##

    ##

    ##

    ##

    ##

    ##

    ##

    ##

    ##

    ##

    ##

    ##
