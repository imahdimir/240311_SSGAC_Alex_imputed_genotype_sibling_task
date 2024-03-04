"""


"""

import pandas as pd
from bgen_reader import open_bgen
from path import Path
from pysnptools.snpreader import Bed
from pysnptools.snpreader import SnpData

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



dyr = Dirs()


def filter_snps() :
    """
    filter out 1000 snps with maf > 0.01 & ( .3 < info <= .31 ) or snps with maf > 0.01 & ( .99 < info <= 1 )
    """
    snps_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/unpacked_ukb_imp_mfi_tgz/ukb_mfi_chr22_v3.txt'

    df = pd.read_csv(snps_fp , sep = '\s+' , header = None)

    # filter snps with maf > 0.01 & ( .3 < info <= .31 )

    # keep snps with maf > 1%
    msk = df[mc.maf].gt(.01)

    # keep snps with ( .3 < info <= .31 )
    msk &= (df[mc.info].gt(.3)) & (df[mc.info].le(.31))

    dfa = df[msk]

    if len(dfa) > 1000 :
        dfa = dfa.sample(1000)

    dfa.to_csv(dyr.snps_2_keep / 'snps_30.txt' ,
               index = False ,
               sep = '\t' ,
               header = False)

    # keep snps with maf > 1%
    msk = df[mc.maf].gt(.01)

    # keep snps with .99 < info <= 1
    msk &= (df[mc.info].gt(.99)) & (df[mc.info].le(1))

    dfa = df[msk]

    if len(dfa) > 1000 :
        dfa = dfa.sample(1000)

    dfa.to_csv(dyr.snps_2_keep + '/snps_99.txt' ,
               index = False ,
               sep = '\t' ,
               header = False)

def filter_sibs():
    """ """

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

    dfb.to_csv(dyr.proj_dir + '/sibs.txt' ,
               index = False ,
               header = False ,
               sep = '\t')

def compute_corr_dosages():
    """"""
    ##
    sibs_fp = "/disk/genetics/ukb/alextisyoung/haplotypes/relatives/bedfiles/hap.kin0"
    df_sibs = pd.read_csv(sibs_fp , sep = '\s+' , dtype = 'string')

    ##
    msk = df_sibs['InfType'].eq('FS')

    df_fs = df_sibs[msk]

    ##
    raw_fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_Tammys_way/bgen_30.raw'
    df_raw = pd.read_csv(raw_fp , sep = '\s+', dtype = 'string')

    ##
    df_raw = df_raw.drop(columns = ['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])

    ##
    df_sib1 = pd.merge(df_fs[['ID1']] , df_raw , left_on = 'ID1' , right_on = 'IID', how = 'left')
    df_sib2 = pd.merge(df_fs[['ID2']] , df_raw , left_on = 'ID2' , right_on = 'IID', how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['IID'])
    df_sib2 = df_sib2.drop(columns = ['IID'])

    ##
    df_sib1 = df_sib1.rename(columns = {'ID1' : 'IID'})
    df_sib2 = df_sib2.rename(columns = {'ID2' : 'IID'})

    ##
    df_sib1 = df_sib1.astype('float')
    df_sib2 = df_sib2.astype('float')

    ##
    df_cors = df_sib1.corrwith(df_sib2, method='pearson')

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
    df_raw = pd.read_csv(raw_fp , sep = '\s+', dtype = 'string')

    ##
    df_raw = df_raw.drop(columns = ['FID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE'])

    ##
    df_sib1 = pd.merge(df_fs[['ID1']] , df_raw , left_on = 'ID1' , right_on = 'IID', how = 'left')
    df_sib2 = pd.merge(df_fs[['ID2']] , df_raw , left_on = 'ID2' , right_on = 'IID', how = 'left')

    ##
    df_sib1 = df_sib1.drop(columns = ['IID'])
    df_sib2 = df_sib2.drop(columns = ['IID'])

    ##
    df_sib1 = df_sib1.rename(columns = {'ID1' : 'IID'})
    df_sib2 = df_sib2.rename(columns = {'ID2' : 'IID'})

    ##
    df_sib1 = df_sib1.astype('float')
    df_sib2 = df_sib2.astype('float')

    ##
    df_cors = df_sib1.corrwith(df_sib2, method='pearson')

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
    from path import Path

    fp = Path(dyr.plink_out / 'bgen_30.bgen')
    fp.exists()

    ##
    fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out/bgen_30.bgen'
    from bgen_reader import open_bgen

    bgen = open_bgen(fp, verbose=False)

    ##
    fp = '/disk/genetics/ukb/mahdimir/imputed_genotype_corr_Tammys_analysis_replication/plink_out_3/bgen_30.bgen'
    bgen = open_bgen(fp, verbose=False)

    ##
    print(bgen.ids)




    ##



    ##



    ##



    ##



    ##



    ##


    ##



    ##
    x = bgen.read(0)

    ##
    x.shape

    ##



    ##
    bgen.samples
    bgen.ids
    bgen.nvariants
    bgen.read(0)

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


        ##
        s1 = df_sib1.iloc[: , [0]]
        s2 = df_sib2.iloc[: , [0]]

        # get pairwise correlation
        s1.corrwith(s2, method='pearson')

        ##
        x = df_sib1.iloc[: , [0]].to_numpy()
        y = df_sib2.iloc[: , [0]].to_numpy()

        ##
        from scipy.stats import pearsonr

        pearsonr(df_sib1.iloc[:, 1] , df_sib2.iloc[:, 1])

        ##
        df_sib1.corrwith(df_sib2, method='pearson')




        ##

        ##
        df = pd.DataFrame(np.random.random((5 , 5)) ,
                       columns = ['gene_' + chr(i + ord('a')) for i in
                                  range(5)])
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
            correlations[col_a + '__' + col_b] = pearsonr(df.loc[:, col_a], df.loc[:, col_b])

        ##








        ##
        df_sib1.corrwith(df_sib2, method='pearson')

        ##
        df_sib1[['rs535733662_G']].corrwith(df_sib2[['rs535733662_G']], method='pearson')

        ##

        df_cors = df_sib1.corrwith(df_sib2, method='pearson')

        ##
        df_sib1.equals(df_sib2)

        ##


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




    ##


    ##


    ##


    ##

    ##





    ##


    ##


    ##


    df = pr.read_bed(bed_fp , as_df = True)

    ##

    df_bed_0 = pd.read_csv(bed_fp , sep= '\t' , comment= 't' , header=None)


    ##

    path = pr.get_example_path("aorta.bed")
    gr = pr.read_bed(path)
    gr

    ##



    ##

    df_bed_0 = pd.read_csv(bed_fp , sep = '\t+' , dtype = 'string')

    ##



    ##

    ##

    ##

def main() :
    pass

    ##
    filter_snps()

    ##

    ##


    ##


    ##
