#!/usr/bin/env python3

import sys
import pandas as pd
import Update_alignments_module as Upalign

def main():

    # Indicate column name given to the projection dataframes
    projection_header = "Projections"

    ## open necessary files

    # open MSA_filenames
    MSA_filenames = Upalign.read_data(sys.argv[1], names=["MSA_filenames"])
    
    # open MSA_transcripts
    MSA_transcripts = Upalign.read_data(sys.argv[2], names=["MSA_transcripts"])
    
    # read included assembly IDs
    assembly_ids = list(Upalign.read_data(sys.argv[3], names=["Assembly_ids"]).loc[:, "Assembly_ids"])
    # No MSA classification file is needed for REFERENCE, so remove REFERENCE from list
    assembly_ids.remove("REFERENCE")

    for assembly in assembly_ids:

        # report progress
        sys.stdout.write(f'\nprocessing {assembly}\n')

        f'{sys.argv[5]}{assembly}_sel_projections_I.txt'

        # open intact projections
        projections_I = Upalign.read_data(f'{sys.argv[4]}{assembly}_sel_projections_I.txt', names=[projection_header])

        # open partial intact projections
        projections_PI = Upalign.read_data(f'{sys.argv[4]}{assembly}_sel_projections_PI.txt', names=[projection_header])

        # open uncertain loss projections
        projections_UL = Upalign.read_data(f'{sys.argv[4]}{assembly}_sel_projections_UL.txt', names=[projection_header])

        # create MSA_classifications file
        MSA_classifications = Upalign.MSA_classification(MSA_filenames, MSA_transcripts, projections_I, 
                                                projections_PI, projections_UL, projection_header=projection_header)

        MSA_classifications.to_csv(f'{sys.argv[5]}{assembly}_MSA_classifications.tsv', sep="\t", index=False)

if __name__ == "__main__":
    main()