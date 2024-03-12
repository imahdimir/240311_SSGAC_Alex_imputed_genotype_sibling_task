"""


    """

import pandas as pd
from bgen_reader import open_bgen
from pathlib import Path
import numpy as np
import itertools

class MfiCols :
    maf = 5
    info = 7

mc = MfiCols()

class Dirs :
    med = '/var/genetics/ws/mahdimir/med/imputed_genotype-sibling-task-240311'
    med = Path(med)
    plink_out = med / 'plink_out'
    
dyr = Dirs()

##
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

    dfa.to_csv(dyr.med / 'snps_30.txt' ,
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

##
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

##
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

##
def make_df_of_iids_from_bgen_open_obj(bg_opn) :
    """ """

    ##
    bg_opn_s = list(bg_opn.samples)
    df = pd.DataFrame({
            'IID' : bg_opn_s
            })

    # get the IID from FID_IID
    df['IID'] = df['IID'].str.split('_').str[1]

    ##
    return df

##
def open_bgen_ret_iid_df_and_prob_arr(bgen_fp) :
    """ """

    ##
    bg_opn = open_bgen(bgen_fp)

    ##
    df_id = make_df_of_iids_from_bgen_open_obj(bg_opn)

    ##
    nd_p = bg_opn.read()

    ##
    return df_id , nd_p

##
def save_dosages_of_all_vars_from_bgen(bgen_fp: Path) :
    """ """

    ##
    df_id , nd_p = open_bgen_ret_iid_df_and_prob_arr(bgen_fp)

    ##
    nd_d = nd_p[: , : , 1] + 2 * nd_p[: , : , 2]

    ##
    df1 = pd.DataFrame(nd_d)

    ##
    df_d = pd.concat([df_id , df1] , axis = 1)

    ##
    _fp = dyr.med / f'dosages_{bgen_fp.stem}.prq'
    df_d.to_parquet(_fp , index = False)

##
def save_hard_calls_of_all_vars_from_bgen(bgen_fp) :
    """ """

    ##
    df_id , nd_p = open_bgen_ret_iid_df_and_prob_arr(bgen_fp)

    ##
    nd_h = np.argmax(nd_p , axis = 2)

    ##
    df1 = pd.DataFrame(nd_h)

    ##
    df_h = pd.concat([df_id , df1] , axis = 1)

    ##
    _fp = dyr.med / f'hard_calls_{bgen_fp.stem}.prq'
    df_h.to_parquet(_fp , index = False)

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

def get_dosages_and_hardcall_of_all_vars_fr_bgen(bgen_fp) :
    """ """

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
    _f = get_dosages_and_hardcall_of_all_vars_fr_bgen
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

    pairs = {
            'sibs'             : '' ,
            'parent_offspring' : '_po' ,
            }

    ##
    prd = itertools.product(info_scores.values() , pairs.values())

    for info , pair in prd :
        print(info , pair)
        fp = dyr.plink_out / f'snps_{info}{pair}.bgen'
        print(fp)

        save_dosages_of_all_vars_from_bgen(fp)

        save_hard_calls_of_all_vars_from_bgen(fp)

    ##

    ##
