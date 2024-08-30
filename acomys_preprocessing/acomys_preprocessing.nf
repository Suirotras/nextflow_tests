params.mountDir = '/home/jari/'

process RETRIEVE_OVERVIEW {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"
    // publishDir 'data/overview_table/', mode: 'copy'

    input:
    // overview.tsv file URL
    val overview_link

    output:
    path "overview.table.tsv"

    """
    #!/usr/bin/env bash

    # retrieve overview.tsv file
    wget ${overview_link} --no-check-certificate -O "overview.table.tsv"
    """
}

process GENERATE_ASSEMBLYLIST {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"
    // publishDir 'data/included_assemblies/', pattern: 'Assembly_names_aBSREL_analysis.csv', mode: 'copy'
    // publishDir 'data/species_df/', pattern: '*_spec_df.csv', mode: 'copy'
    
    input:
    // Incl_assembly_list_generator.R input
    path overview
    
    output:
    // output files
    path 'Assembly_names_aBSREL_analysis.csv', emit: assembly_ids
    path 'dup_spec_df.csv', emit: dup_spec_df
    path 'single_spec_df.csv', emit: single_spec_df

    """
    #!/usr/bin/env bash

    Rscript ${launchDir}/scripts/Incl_assembly_list_generator.R ${overview} 'Assembly_names_aBSREL_analysis.csv' 'dup_spec_df.csv' 'single_spec_df.csv'
    """
}

process RETRIEVE_I_PI_UL_TRANSCRIPTS {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"

    input:
    // base URL
    val baselink
    // list of assembly IDs
    path assembly_ids
    // codon alignments URL
    val multiple_alignment_link
    // gene loss summary file
    val loss_summ_link

    output:
    path 'intermediate_files/MSA_transcripts.txt', emit: MSA_transcripts
    path 'intermediate_files/MSA_filenames.txt', emit: MSA_filenames
    path 'sel_transcripts_projections', emit: sel_transcripts

    """
    #!/usr/bin/env bash

    mkdir -p intermediate_files
    mkdir -p sel_transcripts_projections

    python3 ${launchDir}/scripts/retrieve_I_PI_UL_transcripts.py ${baselink} ${assembly_ids} ${multiple_alignment_link} ${loss_summ_link}
    """
}

process GET_MSA_CLASSIFICATIONS {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"

    input:
    path MSA_filenames
    path MSA_transcripts
    path assembly_ids
    path sel_transcripts

    output:
    path 'MSA_classifications', emit: MSA_class

    """
    #!/usr/bin/env bash

    mkdir -p MSA_classifications

    python3 ${launchDir}/scripts/Get_MSA_classifications.py ${MSA_filenames} ${MSA_transcripts} ${assembly_ids} ${sel_transcripts}/ "MSA_classifications/"
    """
}

process GET_PROJECTION_FASTA {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"

    input:
    path MSA_class
    val acomys_alignment_link
    path fully_intact

    output:
    path 'acomys_seq_fasta', emit: acomys_seq_fasta
    path 'acomys_added_projections/HLacoCah2_added_projections.tsv', emit: acomys_added_proj

    """
    #!/usr/bin/env bash

    mkdir -p ./acomys_seq_fasta
    mkdir -p ./acomys_added_projections

    python3 ${launchDir}/scripts/Get_projection_fasta.py ${MSA_class}/HLacoCah2_MSA_classifications.tsv ${acomys_alignment_link} ${fully_intact} "acomys_seq_fasta/" "acomys_added_projections/HLacoCah2_added_projections.tsv"
    """
}

// Temporary process to limit the number of alignments that are processed
process REDUCE_ADDED_PROJECTIONS {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"

    input:
    path acomys_added_proj

    output:
    path 'Added_projections_small_sample.tsv', emit: small_sample

    """
    #!/usr/bin/env bash

    head -n 20 ${acomys_added_proj} > 'Added_projections_small_sample.tsv' 
    """
}

process MACSEV2_ENRICHALIGNMENT {
    container 'jarivdiermen/java_gnuparallel_ubuntu:latest'
    containerOptions "-B $params.mountDir:$params.mountDir"

    input:
    path added_acomys_projections
    val multiple_alignment_link
    path fasta_output
    val MACSEv2_URL

    output:
    path 'cactus.alignments_updated_stats', emit: updated_stats
    path 'cactus.alignments_updated_NT', emit: updated_NT
    path 'cactus.alignments_updated_AA', emit: updated_AA

    """
    #!/usr/bin/env bash

    mkdir -p cactus.alignments
    mkdir -p cactus.alignments_updated_stats
    mkdir -p cactus.alignments_updated_NT
    mkdir -p cactus.alignments_updated_AA

    bash ${launchDir}/scripts/MACSEv2_enrichAlignment.sh ${added_acomys_projections} ${multiple_alignment_link} cactus.alignments/ ${fasta_output}/ ${MACSEv2_URL}

    # move alignment files to new folders
    mv cactus.alignments/*_stats.csv "cactus.alignments_updated_stats"
    mv cactus.alignments/*_NT.fasta "cactus.alignments_updated_NT"
    mv cactus.alignments/*_AA.fasta "cactus.alignments_updated_AA"
    """
}


baselink = channel.value("https://genome.senckenberg.de/download/TOGA/human_hg38_reference/")
// URL for overview.tsv file
overview_link = channel.value("https://genome.senckenberg.de/download/TOGA/human_hg38_reference/overview.table.tsv")
// URL for multipleCodonAligment files
multiple_alignment_link = channel.value("https://genome.senckenberg.de/download/TOGA/human_hg38_reference/MultipleCodonAlignments/")
// URL for codon alignments (i.e. updated CDS against hg38 reference CDS)
acomys_alignment_link = channel.value("https://genome.senckenberg.de/download/TOGA/human_hg38_reference/Rodentia/Acomys_cahirinus__Egyptian_spiny_mouse__HLacoCah2/codonAlignments.fa.gz")
// Added to create the gene loss summary URL
loss_summ_link="loss_summ_data.tsv.gz"

// Path to bed file with fully intact HLacoCah2 transcripts
aco_fully_intact = channel.value("${launchDir}/input_files/HLacoCah2.query_annotation.fullyIntact.bed.gz")

// URL to the MACSEv2 program
MACSEv2_URL = channel.value("https://www.agap-ge2pop.org/wp-content/uploads/macse/releases/macse_v2.07.jar")

workflow {

    // get overview.tsv file
    overview_ch = RETRIEVE_OVERVIEW(overview_link)
    //overview_ch.flatten().view()

    // Generate the genome assembly list
    assemblylist_ch = GENERATE_ASSEMBLYLIST(overview_ch.flatten())

    // assemblylist_ch.assembly_ids.view()
    // assemblylist_ch.dup_spec_df.view()
    // assemblylist_ch.single_spec_df.view()

    // Generate the transcript files
    transcripts_ch = RETRIEVE_I_PI_UL_TRANSCRIPTS(baselink, assemblylist_ch.assembly_ids, multiple_alignment_link, loss_summ_link)

    //transcripts_ch.MSA_transcripts.view()
    //transcripts_ch.MSA_filenames.view()
    //transcripts_ch.sel_transcripts.view()

    MSA_class_ch = GET_MSA_CLASSIFICATIONS(transcripts_ch.MSA_filenames, transcripts_ch.MSA_transcripts, assemblylist_ch.assembly_ids, transcripts_ch.sel_transcripts)
    //MSA_class_ch.MSA_class.view()

    get_proj_ch = GET_PROJECTION_FASTA(MSA_class_ch.MSA_class, acomys_alignment_link, aco_fully_intact)

    reduced_ch = REDUCE_ADDED_PROJECTIONS(get_proj_ch.acomys_added_proj)

    macsev2_ch = MACSEV2_ENRICHALIGNMENT(reduced_ch.small_sample, multiple_alignment_link, get_proj_ch.acomys_seq_fasta, MACSEv2_URL)

}