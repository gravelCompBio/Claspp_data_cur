import getData.DataProcess as dp



text_loc="data/text_loc/"


def main():
    lab2peps, lab2affm, lab2maxnc = dp.prep_data_for_spec_clus(text_loc=text_loc)
    dp.spec_clus_sk_learn(lab2peps,lab2affm,lab2maxnc, text_loc)





if __name__ == "__main__":
    main()