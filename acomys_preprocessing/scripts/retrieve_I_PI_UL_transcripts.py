#!/usr/bin/env python3

import sys
import pandas as pd
import Update_alignments_module as Upalign

def main():

    ## retrieve list of URLs where TOGA data of every species is stored
    baselink = sys.argv[1]
    # Retrieve URLs from baselink. excl_filters and incl_filters remove URLs where unneeded data is stored.
    #baselink_URLs = Upalign.retrieve_options(baselink, excl_filters=["?", "Matrix/", "MultipleCodonAlignments/", 
    #                                                     	"hg38_chains/", "overview.table.tsv"], 
    #                                            incl_filters=["human_hg38_reference/Afrotheria"])
    baselink_URLs = Upalign.retrieve_options(baselink, excl_filters=["?", "Matrix/", "MultipleCodonAlignments/", 
                                                        "hg38_chains/", "overview.table.tsv"], 
                                                incl_filters=["human_hg38_reference"])
    
    # List to store all URLs
    all_URLs = []
    # Process subsequent URLs recursively
    for URL in baselink_URLs:
        sub_URLs = Upalign.retrieve_options(URL, excl_filters=["?"], incl_filters=baselink_URLs)
        all_URLs.extend(sub_URLs)
    
    # read included assembly IDs
    assembly_ids = list(Upalign.read_data(sys.argv[2], names=["Assembly_ids"]).loc[:, "Assembly_ids"])
    
    # Only keep URLs that store data on included species
    included_URLs = Upalign.filter_URLs_by_substrings(all_URLs, assembly_ids)

    # Retrieve reference transcripts used for MSAs
    MSA_filenames, MSA_transcripts = Upalign.MSA_filename_transcript_retrieval_url(sys.argv[3])

    # retrieve selected projections and transcripts for all included species
    for URL in included_URLs:
        # report progress
        sys.stdout.write(f'\nprocessing {URL}\n')

        # retrieve gene loss summary file
        loss_summary = Upalign.gene_loss_summary_retrieval(URL + sys.argv[4])
    
        # remove gene IDs from loss summary file
        loss_summary = loss_summary[~loss_summary.loc[:, "Type"].isin(["GENE"])]

        # Filter for MSA_transcripts
        MSA_loss_sum = Upalign.filter_summary_file_transcripts(loss_summary, MSA_transcripts)

        # Filter for intact, partial intact AND uncertain loss transcripts
        MSA_sum_sel = Upalign.filter_summary_file_classification(MSA_loss_sum, classification = ["I", "PI", "UL"])
        # filter for projections
        sel_projections = MSA_sum_sel[MSA_sum_sel.loc[:, "Type"] == "PROJECTION"].loc[:, "Gene_Transcript"]
        # filter for transcripts
        sel_transcripts = MSA_sum_sel[MSA_sum_sel.loc[:, "Type"] == "TRANSCRIPT"].loc[:, "Gene_Transcript"]
    
        # for loop creating dictionary with Intact, partial intact or uncertain loss projections and transcripts
        sel_projections_dict = {}
        sel_projections_dict["I_PI_UL"] = sel_projections
        sel_transcripts_dict = {}
        sel_transcripts_dict["I_PI_UL"] = sel_transcripts
    
        for classification in ["I", "PI", "UL"]:
            MSA_sum_sel_filtered = Upalign.filter_summary_file_classification(MSA_loss_sum, classification = [classification])
            # filter for projections and add to sel_projections_dict
            sel_projections_filtered = MSA_sum_sel_filtered[MSA_sum_sel_filtered.loc[:, "Type"] == "PROJECTION"].loc[:, "Gene_Transcript"]
            sel_projections_dict[classification] = sel_projections_filtered
            # filter for transcripts and add to sel_transcripts_dict
            sel_transcripts_filtered = MSA_sum_sel_filtered[MSA_sum_sel_filtered.loc[:, "Type"] == "TRANSCRIPT"].loc[:, "Gene_Transcript"]
            sel_transcripts_dict[classification] = sel_transcripts_filtered
    
        URL_assembly_id = URL.split("__")[2][:-1]

        # write the projections from dictionary to their respective files
        for key, value in sel_projections_dict.items():
            with open(f'./sel_transcripts_projections/{URL_assembly_id}_sel_projections_{key!s}.txt', mode="w", encoding="utf-8") as f:
                Upalign.write_series(f, value)
    
        # write the transcripts from dictionary to their respective files
        for key, value in sel_projections_dict.items():
            with open(f'./sel_transcripts_projections/{URL_assembly_id}_sel_transcripts_{key!s}.txt', mode="w", encoding="utf-8") as f:
                Upalign.write_series(f, value)

    # write the MSA_transcripts to file
    with open("./intermediate_files/MSA_transcripts.txt", mode="w", encoding="utf-8") as f:
        Upalign.write_series(f, MSA_transcripts)

    # write the MSA_filenames to file
    with open("./intermediate_files/MSA_filenames.txt", mode="w", encoding="utf-8") as f:
        Upalign.write_series(f, MSA_filenames)

if __name__ == "__main__":
    main()