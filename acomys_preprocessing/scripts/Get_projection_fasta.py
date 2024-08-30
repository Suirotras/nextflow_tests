#!/usr/bin/env python3

import sys
import pandas as pd
import pyfastx
import Update_alignments_module as Upalign

def main():
    
    # define codon_alignment.fa as name for the downloaded fasta file
    name = "codonAlignments.fa.gz"

    #First download the fasta file from server
    Upalign.download_gzipped_file(sys.argv[2], name=name)

    # retrieve fully intact HLacoCah2 projections
    fully_intact = list(pd.read_table(sys.argv[3], compression="gzip").iloc[:, 3])
    
    # open fasta file
    fa = pyfastx.Fasta(name, full_name = True)
    # open MSA_classifications
    MSA_classifications = pd.read_table(sys.argv[1])
    # select MSA_transcript column
    MSA_transcripts = MSA_classifications.loc[:, "MSA_transcript"]
    # remove duplicate transcripts
    MSA_transcripts = MSA_transcripts.drop_duplicates()

    # list of projections selected for addition to the MSAs
    added_projections_list = []
    
    # loop through MSA_transcripts
    for transcript in MSA_transcripts:
        # report progress
        sys.stdout.write(f'\nprocessing {transcript}\n')

        # filter for MSA_transcript
        MSA_class_filtered = MSA_classifications[MSA_classifications.loc[:, "MSA_transcript"] == transcript]
        # find projections for intact, partial intact or uncertain loss.
        Projection_dict = {}
        for classification in ["I", "PI", "UL"]:
            Projection_dict[classification] = MSA_class_filtered[MSA_class_filtered.loc[:, "Classification"] == classification].loc[:, "projection"]

        for classification, projection in Projection_dict.items():
            
            n_of_projections = len(projection)

            # if zero projections have this classification, skip to next classification
            if n_of_projections == 0:
                continue

            if classification == "I":
                # if one intact projection, write this to fasta
                if n_of_projections == 1:
                    # retrieve the chosen projection 
                    seq = Upalign.select_sequence(fa, fa.keys(), identifier=list(projection)[0], select="QUERY")
                    
                    # write fasta
                    Upalign.write_seq_to_fasta(seq, file=f'./{sys.argv[4]}HLacoCah2_{transcript}.fa',
                                       seq_header=f'>vs_HLacoCah2\t{list(projection)[0]}')
                    
                    #print info to terminal
                    sys.stdout.write(f'One "{classification}" projection: {list(projection)[0]}\n')
                    added_projections_list.append(list(projection)[0])
                    break 
                    
                if n_of_projections > 1:
                    # check how many are fully intact
                    fully_intact_projections = projection[projection.isin(fully_intact)]
                    # if one fully intact left, write to fasta
                    if len(fully_intact_projections) == 1:
                        # retrieve the chosen projection 
                        seq = Upalign.select_sequence(fa, fa.keys(), identifier=list(fully_intact_projections)[0], select="QUERY")

                        # write fasta
                        Upalign.write_seq_to_fasta(seq, file=f'./{sys.argv[4]}HLacoCah2_{transcript}.fa',
                                       seq_header=f'>vs_HLacoCah2\t{list(fully_intact_projections)[0]}')
                        
                        #print info to terminal
                        sys.stdout.write(f'Multiple "{classification}" projections: {list(projection)}\nOne fully "{classification}" projection: {list(fully_intact_projections)[0]}\n')
                        added_projections_list.append(list(fully_intact_projections)[0])
                        break

                    elif len(fully_intact_projections) > 1:
                        # retrieve projection similarity scores between QUERY and REFERENCE
                        similarity_scores = Upalign.QUERY_REFERENCE_similarity_scores(fa, fully_intact_projections)
                        # get projection with highest similarity score
                        max_similarity_projection = max(similarity_scores, key=similarity_scores.get)
                        # retrieve and write chosen projection
                        seq = Upalign.select_sequence(fa, fa.keys(), identifier=max_similarity_projection, select="QUERY")
                        Upalign.write_seq_to_fasta(seq, file=f'./{sys.argv[4]}HLacoCah2_{transcript}.fa',
                                       seq_header=f'>vs_HLacoCah2\t{max_similarity_projection}')
                        
                        #print info to terminal
                        sys.stdout.write(f'Multiple "{classification}" projections: {list(projection)}\nMultiple fully "{classification}" projections: {list(fully_intact_projections)}\nMax similarity projection: {max_similarity_projection}\n')
                        added_projections_list.append(max_similarity_projection)
                        break
                                
            if classification in ["PI", "UL"]:
                # if one partial intact or uncertain loss projection, write to fasta
                if n_of_projections == 1:
                    # retrieve the chosen projection 
                    seq = Upalign.select_sequence(fa, fa.keys(), identifier=list(projection)[0], select="QUERY")
                    
                    # write fasta
                    Upalign.write_seq_to_fasta(seq, file=f'./{sys.argv[4]}HLacoCah2_{transcript}.fa',
                                       seq_header=f'>vs_HLacoCah2\t{list(projection)[0]}')
                    
                    #print info to terminal
                    sys.stdout.write(f'One "{classification}" projection: {list(projection)[0]}')    
                    added_projections_list.append(list(projection)[0])
                    break
                
                if n_of_projections > 1:
                    # retrieve projection similarity scores between QUERY and REFERENCE
                    similarity_scores = Upalign.QUERY_REFERENCE_similarity_scores(fa, projection)
                    # get projection with highest similarity score
                    max_similarity_projection = max(similarity_scores, key=similarity_scores.get)
                    # retrieve and write chosen projection
                    seq = Upalign.select_sequence(fa, fa.keys(), identifier=max_similarity_projection, select="QUERY")
                    Upalign.write_seq_to_fasta(seq, file=f'./{sys.argv[4]}HLacoCah2_{transcript}.fa',
                                    seq_header=f'>vs_HLacoCah2\t{max_similarity_projection}')
                    
                    #print info to terminal
                    sys.stdout.write(f'Multiple "{classification}" projections: {list(projection)}\nMax similarity projection: {max_similarity_projection}\n') 
                    added_projections_list.append(max_similarity_projection)
                    break
    
    # Record the projections selected for MSA addition to file
    added_projections = MSA_classifications[MSA_classifications.loc[:, "projection"].isin(added_projections_list)]
    added_projections.to_csv(sys.argv[5], sep="\t", index=False)

if __name__ == "__main__":
    main()