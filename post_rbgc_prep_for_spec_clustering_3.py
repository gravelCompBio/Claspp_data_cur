import getData.DataProcess as dp



text_loc="data/text_loc/"


def main():
    masterlist, pep2clus = dp.readin_files_v1(text_loc = text_loc)
    pep2info, hpep2labs, bs, mlabpeps, pep2mllabel, prank = dp.track_all_human_multi_lable_peps(masterlist=masterlist,text_loc=text_loc)
    masterlist, nonIntrestClusters, clus2info, nonHumanClus, rank, mlinfokey = dp.set_sampeling_prority(masterlist=masterlist, mlabpeps=mlabpeps, pep2mllabel=pep2mllabel, prank=prank, pep2info=pep2info, bs=bs)
    dp.set_up_medium_negitive_data_set(nonIntrestClusters=nonIntrestClusters, clus2info=clus2info, nonHumanClus=nonHumanClus, rank=rank , text_loc=text_loc)
    dp.set_up_for_spec_clustering(masterlist=masterlist,mlinfokey=mlinfokey,text_loc=text_loc)
    dp.prep_matlab_dist_mat_for_USSC(text_loc=text_loc)

if __name__ == "__main__":
    main()