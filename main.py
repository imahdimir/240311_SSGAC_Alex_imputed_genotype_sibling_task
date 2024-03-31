"""


    """

import itertools
import numpy as np
import pandas as pd
from bgen_reader import open_bgen
from pathlib import Path

class MfiCols :
    maf = 5
    info = 7

mc = MfiCols()

class Dirs :
    med = '/var/genetics/ws/mahdimir/med/imputed_genotype-sibling-task-240311'
    med = Path(med)
    plink_out = med / 'plink_out'
    out = '/var/genetics/ws/mahdimir/DropBox/0-all/1-out-all/imputed_genotype-sibling-task-240311'
    out = Path(out)
    out_dta = out / 'dta'

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

    dfb.to_csv(dyr.med + '/sibs.txt' ,
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

    # keep only parent offspring pairs
    msk = df['InfType'].eq('PO')
    dfa = df[msk]

    df1 = dfa[['FID1' , 'ID1']]
    df2 = dfa[['FID2' , 'ID2']]

    df1.columns = ['FID' , 'ID']
    df2.columns = ['FID' , 'ID']

    dfb = pd.concat([df1 , df2])
    dfb = dfb.iloc[: , :2]

    dfb.to_csv(dyr.med + '/parent_offspring.txt' ,
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
def gat_pairs_ids_in_pairs(identifier) :
    """ """

    rel_fp = '/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0'
    df = pd.read_csv(rel_fp , sep = '\s+' , dtype = 'string')

    msk = df['InfType'].eq(identifier)

    df = df[msk]

    df = df[['ID1' , 'ID2']]

    return df

##
def make_prq_fp(gts_type , info_score , pair_suf) :
    """ """
    prq_fp = dyr.med / f'{gts_type}_snps_{info_score}{pair_suf}.prq'
    return prq_fp

##
def make_pairs_gts_dfs(df_gts , df_pairs_ids) :
    """ """

    ##
    dfa = pd.merge(df_pairs_ids[['ID1']] ,
                   df_gts ,
                   left_on = 'ID1' ,
                   right_on = 'IID' ,
                   how = 'left')
    dfb = pd.merge(df_pairs_ids[['ID2']] ,
                   df_gts ,
                   left_on = 'ID2' ,
                   right_on = 'IID' ,
                   how = 'left')

    ##
    dfa = dfa.drop(columns = ['ID1'])
    dfb = dfb.drop(columns = ['ID2'])

    ##
    return dfa , dfb

##
def save_corr_of_sib_pairs_with_gts_df(pair_identifier ,
                                       gts_type ,
                                       info_score ,
                                       pair_type ,
                                       pair_suf
                                       ) :
    """ """

    ##
    df_pairs_ids = gat_pairs_ids_in_pairs(pair_identifier)

    ##
    prq_fp = make_prq_fp(gts_type , info_score , pair_suf)
    df_gts = pd.read_parquet(prq_fp)

    ##
    dfa , dfb = make_pairs_gts_dfs(df_gts , df_pairs_ids)

    ##
    gts1 = dfa.iloc[: , 1 :]
    gts2 = dfb.iloc[: , 1 :]

    ##
    df_cors = gts1.corrwith(gts2 , method = 'pearson')

    ##
    out_fp = dyr.out_dta / f'corr_{pair_type}_{gts_type}_{info_score}.xlsx'
    df_cors.to_excel(out_fp , index = False , header = False)
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
            'sibs'             : ('' , 'FS') ,
            'parent_offspring' : ('_po' , 'PO') ,
            }

    gts_types = {
            0 : 'dosages' ,
            1 : 'hard_calls'
            }

    ##
    prd = itertools.product(info_scores.values() , pairs.values())

    for info , pair in prd :
        print(info , pair)
        fp = dyr.plink_out / f'snps_{info}{pair[0]}.bgen'
        print(fp)

        save_dosages_of_all_vars_from_bgen(fp)

        save_hard_calls_of_all_vars_from_bgen(fp)

    ##
    prd = itertools.product(info_scores.values() ,
                            pairs.keys() ,
                            gts_types.values())

    for info , pair_type , gm in prd :
        pair_iden = pairs[pair_type][1]
        print(info , pair_type , gm)
        pair_suf = pairs[pair_type][0]
        save_corr_of_sib_pairs_with_gts_df(pair_iden ,
                                           gm ,
                                           info ,
                                           pair_type ,
                                           pair_suf
                                           )


    ##




    ##

    ##

    ##

    ##
