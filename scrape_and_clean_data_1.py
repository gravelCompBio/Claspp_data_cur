import os
import time
import getData.WebSrcape as ws
import getData.DataProcess as dp
import getData.FileReader as fr



dbdir = "data/dbPTMloc/"
csv_loc = "data/csv_loc/"
fasta_loc = "data/fasta_loc/"  
text_loc = "data/text_loc/"





def main():
    if not os.path.isdir(f"{dbdir}"):
        os.makedirs(f"{dbdir}")
    if not os.path.isdir(f"{fasta_loc}"):
        os.makedirs(f"{fasta_loc}")
    
    print("step 1) downloading the following PTM types")
    ws.download_db_PTM()
    
    time.sleep(15)
    

    print("step 2) reading and intial cleaning")
    uni2name , masterlist = dp.unzip_and_read_files_from_dbPTM(db_loc=dbdir)
    if len(masterlist)==0:
        "weird run again"
        return 
    missingpeps, missingunis, uni2seq, masterlist = ws.find_missing_data(masterlist=masterlist)
    masterlist, uniset = dp.intial_hand_curation_update(masterlist=masterlist, csv_loc=csv_loc)

    print("step 3) fining missing data points in dbPTM according to uniprot database standards")
    if not os.path.exists(f"{fasta_loc}uniprotPostMissingPepPreConRes.fasta"):
        uni2seq = ws.download_fasta_from_unis_v1( uniset=uniset, fasta_loc=fasta_loc)
    else:
        uni2seq = fr.read_downloaded_fasta_from_unis_v1(fasta_loc=fasta_loc)
    newunis, oldunis = dp.separate_old_and_new_data_v1(uni2seq=uni2seq, masterlist=masterlist, text_loc=text_loc)
    print("step 4) mapping and out-of-date data with older versions of uniprot (round 1 of interation)")

    if not os.path.exists(f"{csv_loc}uni2uniparc.csv"):
        uni2uniparc = ws.map_out_uniprotIDs_to_uniParc_v1(oldunis=oldunis, csv_loc=csv_loc)
    else:
        uni2uniparc = fr.read_uniprotIDs_to_uniParc_v1(csv_loc=csv_loc)

    print("step 5) retriving out-of-date sequences (round 1 of interation)")
    if not os.path.exists(f"{fasta_loc}oldunis.fasta"):
        uniparc2uni,olduni2seq = ws.download_fasta_from_parc_v1(uni2uniparc=uni2uniparc, fasta_loc=fasta_loc)
    else:
        uniparc2uni,olduni2seq = fr.read_fasta_from_parc_v1(uni2uniparc=uni2uniparc, fasta_loc=fasta_loc)

    print("step 6) mapping and out-of-date data with older versions of uniprot (round 2 of interation)")
    if not os.path.exists(f"{csv_loc}uni2uniparcV2.csv"):
        uni2uniparc = ws.map_out_uniprotIDs_to_uniParc_v2(oldunis=oldunis, csv_loc=csv_loc)
    else:
        uni2uniparc = fr.read_uniprotIDs_to_uniParc_v2(csv_loc=csv_loc)
    newunis, oldunis = dp.separate_old_and_new_data_v2(uni2seq=uni2seq,olduni2seq=olduni2seq,masterlist=masterlist)

    print("step 7) retriving out-of-date sequences (round 2 of interation)")
    if not os.path.exists(f"{fasta_loc}oldunisV2.fasta"):
        uniparc2uni, olduni2seq = ws.download_fasta_from_parc_v2( uni2uniparc=uni2uniparc, fasta_loc=fasta_loc)
    else:
        uniparc2uni, olduni2seq = fr.read_downloaded_fasta_from_unis_v2(uni2uniparc=uni2uniparc ,fasta_loc=fasta_loc)

    print("step 8) scanning for broken data points found in dbPTM")
    uni2seq = dp.final_merge(uni2seq=uni2seq, olduni2seq=olduni2seq)
    masterlist = dp.remove_broken_peptides(masterlist=masterlist)
    missmatchedseq, updatedlist = dp.fix_missmatched_uniprot_postions(uni2seq=uni2seq, masterlist=masterlist)

    print("step 9) finding out-of-date sequences that change due to fixing")
    if not os.path.exists(f"{fasta_loc}mmUniFull.fasta"):
        ws.download_fasta_from_parc_v3_multiThread(missmatchedseq=missmatchedseq, fasta_loc=fasta_loc)
    mmUni2ver2seq = fr.read_downloaded_fasta_from_unis_v3(fasta_loc=fasta_loc)
    masterlist = dp.fix_missmatched_uniprot_postions_v2( uni2seq=uni2seq, mmUni2ver2seq=mmUni2ver2seq, missmatchedseq=missmatchedseq, updatedlist=updatedlist)
    dp.filter_for_human_and_major_res_per_PTM_type(masterlist=masterlist, text_loc=text_loc)
    print("Finished Cleaning (run the next python file !!!!)")

    






if __name__ == "__main__":
    main()
