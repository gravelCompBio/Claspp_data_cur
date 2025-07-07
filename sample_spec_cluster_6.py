
import os
import getData.DataProcess as dp

text_loc = "data/text_loc/"
csv_loc = "data/csv_loc/"


def main():
    pepsused, mLcluslab2pd, cluslab2pd, cluslabs, cluslab2peps, notclustered = dp.sample_from_spec_clus(csv_loc = csv_loc, text_loc = text_loc)
    negres2peps, negSet = dp.sample_negtive_data(csv_loc=csv_loc, text_loc=text_loc)
    dp.data_split_train_test_val(pepsused=pepsused, mLcluslab2pd=mLcluslab2pd, cluslab2pd=cluslab2pd, negres2peps=negres2peps, negSet=negSet, cluslabs=cluslabs, cluslab2peps=cluslab2peps, notclustered=notclustered, csv_loc=csv_loc)
    dp.reduce_neg_ratio(csv_loc)



if __name__ == "__main__":
    main()