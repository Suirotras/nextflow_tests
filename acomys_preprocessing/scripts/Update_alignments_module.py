# Import modules
import sys
import os
import re
import gzip
import requests
import io
import pandas as pd
from bs4 import BeautifulSoup
from urllib.parse import urljoin
from urllib3.exceptions import InsecureRequestWarning
from urllib3 import disable_warnings
from Bio import Align

def change_working_directory(path):
    """Function to change working directory.

    args:
        path: the path to the new working directory.
    """
    os.chdir(path)

def retrieve_options(url, excl_filters=None, incl_filters=None, verify=False):
    """Function to retrieve list of all sub-URLs at a given URL 
    (e.g. retrieve ["https://hello.world.com/hello/", "https://hello.world.com/world/"] from
    "https://hello.world.com/")

    args:
        url: URL to search
        excl_filters: list of strings. Removes all URLs containing at least one of the strings.
        incl_filters: list of strings. Only keep URLs containing at least one of the strings.
        verify: a boolean. If False, verification of GET request is turned off.
    """
    # Send a GET request to the URL
    disable_warnings(InsecureRequestWarning)
    response = requests.get(url, verify=verify)
    # Create a BeautifulSoup object with the response content
    soup = BeautifulSoup(response.content, "html.parser")
    # Find all <a> tags on the page
    links = soup.find_all("a")
    # Extract the href attribute from each <a> tag
    options = []
    for link in links:
        href = link.get("href")
        if href:
            absolute_url = urljoin(url, href)
            # if exclusion and inclusion criteria selected
            if excl_filters and incl_filters:
                if not any(filter_option in absolute_url for filter_option in excl_filters):
                    if any(filter_option in absolute_url for filter_option in incl_filters):
                        options.append(absolute_url)
            # if only excl_filters selected
            elif excl_filters:
                if not any(filter_option in absolute_url for filter_option in excl_filters):
                    options.append(absolute_url)
            # if only incl_filters selected
            elif incl_filters:
                if any(filter_option in absolute_url for filter_option in incl_filters):
                    options.append(absolute_url)
            # if no filters selected
            else:
                options.append(absolute_url)
    return options

def filter_URLs_by_substrings(string_list, substring_list):
    """Function that takes a list of strings, and returns a list containing only the strings where 
    at least one of the substrings was present.

    args:
        string_list: The list of strings to be filtered
        substring_list: The list of substrings used to filter the list of strings. Only strings containing
        at least one substring are returned.
    """
    filtered_strings = [string for string in string_list for substring in substring_list if re.search(substring, string)]
    return filtered_strings

def MSA_filename_transcript_retrieval_url(MSA_URL):
    """Function to retrieve MSA filename and corresponding transcript IDs.
    Returns two lists, like ["MSA_filename_1", "MSA_filename_2", ...] and 
    ["MSA_transcript_1", "MSA_transcript_ID_2", ...]

    args:
        MSA_URL: URL containing all the multiple sequence alignments (MSAs).
    """
    # disable verification warning
    disable_warnings(InsecureRequestWarning)
    # retrieve page
    MSA_page = web_response = requests.get(MSA_URL, verify=False)
    # get HTML code
    MSA_HTML = web_response.text
    # retrieve filenames and transcript IDs from HTML code
    MSA_regex = re.compile(r'href="((ENST.+?)\..+?\..+?\.gz)"')
    MSA_tuple_list = MSA_regex.findall(MSA_HTML)
    # unpack list of tuples into two lists
    MSA_filenames, MSA_transcripts = zip(*MSA_tuple_list)

    return MSA_filenames, MSA_transcripts

def gene_loss_summary_retrieval(loss_summ_URL):
    """Function to retrieve gene_loss_summary_file from url

    args:
        loss_summ_URL: URL containing the loss_summ_data.tsv.gz file.
        The URL already containts the loss_summary file name.
        (e.g. https://URL_to_summary_file/loss_summ_data.tsv.gz)
    """
    # disable verification warning
    disable_warnings(InsecureRequestWarning)
    # request summary file
    web_response = requests.get(loss_summ_URL, verify=False)
    # retrieve contents from web_response
    loss_summ_gz_file = web_response.content
    # convert to io.BytesIO class object
    loss_summ_gz_file = io.BytesIO(loss_summ_gz_file)
    # Unzip and read summary file
    with gzip.GzipFile(fileobj=loss_summ_gz_file) as f:
        # Passing a binary file to csv.reader works in PY2
        loss_summ_file = pd.read_table(f, header=None, names=["Type", "Gene_Transcript", "Classification"])

    return loss_summ_file

def filter_summary_file_transcripts(loss_sum, transcripts):
    """Function to filter loss_summ_data.tsv file with a list (or pd.series) of transcript IDs.

    args:
        loss_sum: pd.Dataframe of loss_summ_data.tsv file.
        transcripts: list of transcript IDs for filtering. Only projections and transcripts with
        transcript IDs present in this list are returned.
    """
    # split column to allow filtering using transcript IDs only. 
    # (i.e. "ENST00000624481.CD99" becomes "ENST00000624481" and "CD99")
    All_transcripts = loss_sum["Gene_Transcript"].str.split(".", n=1).str[0]
    filtered_loss_sum = loss_sum[All_transcripts.isin(transcripts)]
    # return gene loss summary file containing only the given transcript IDs
    return filtered_loss_sum

def filter_summary_file_classification(loss_sum, classification = ["I", "PI", "UL"]):
    """Function to filter loss_summ_data.tsv file for gene loss classifications.

    args:
        loss_sum: pd.Dataframe of loss_summ_data.tsv file.
        classification: list of gene loss classifications. Only Projections, Transcripts
        and Genes with a classification present in this list are returned. 
    """
    loss_sum_selected = loss_sum[loss_sum.loc[:, "Classification"].isin(classification)]
    return loss_sum_selected

def write_series(file, series):
    """Function to write a list (or pandas.series) to a file.

    args:
        file: the filename to be written, including its file path.
        series: the list, or pandas.series, to be written to the file.
    """
    for row in series:
        file.write(row + "\n")

def read_data(file, names, seperator = ",", header = None):
    """Read in a file that contains columns of data. Return this data
    as a pandas.Dataframe.

    args:
        file: the filename, including its file path, to be read.
        names: The names to give the columns. Should be the same length
        as the number of columns.
        seperator: The separator used to distinguish the columns,
        Defaults to a comma (",").
        header: Indicates whether there is a header row with column names. 
        Defaults to None, indicating no header row.
    """
    return pd.read_csv(file, sep=seperator, header = header, names = names)

def MSA_classification(MSA_filenames, MSA_transcripts, projections_I, projections_PI, projections_UL,
                       projection_header = "Projections"):
    """Function that creates MSA_classification files. These are files indicating for every MSA transcript the following:
        - The projections associated with that MSA transcript.
        - The classification of that projection (Intact, partial intact or uncertain loss).

    Input: takes in the MSA_filenames, MSA_transcripts, and the intact, partial intact and uncertain loss 
    projection pd.dataframes. projection_header represents the column header given to the projection dataframes.
    Output: tab-delimited MSA_classification file relating MSA_transcripts to their corresponding projections. Also lists the 
    classification for each projection (intact (I), partial intact (PI), uncertain loss(UL)).
    """
    # Add transcript ID column by splitting projection IDs
    for projection_data in [projections_I, projections_PI, projections_UL]:
        projection_data["Transcript_ID"] = projection_data.loc[:, projection_header].str.split(".", n=1).str[0]

    MSA_filename_list = []
    MSA_transcript_list = []
    Projection_list = []
    Classification_list = []
    
    for MSA_filename, MSA_transcript in zip(list(MSA_filenames.iloc[:, 0]), list(MSA_transcripts.iloc[:, 0])):
        # for every MSA_filename and transcript, look for projections in projections_I, projections_PI
        # and projections_UL
        for classification, projection_data in zip(["I", "PI", "UL"], [projections_I, projections_PI, projections_UL]):
            # retrieve projections for MSA transcript
            MSA_projections = projection_data[projection_data.loc[:, "Transcript_ID"] == MSA_transcript].loc[:, projection_header]
            # check if projections were found
            if len(MSA_projections) != 0:
                # if one or more projections found, add to list
                for projection in MSA_projections:
                    MSA_filename_list.append(MSA_filename)
                    MSA_transcript_list.append(MSA_transcript)
                    Projection_list.append(projection)
                    Classification_list.append(classification)

    # create MSA_classifications file
    MSA_classifications = pd.DataFrame(data={"MSA_filename": MSA_filename_list, 
                                             "MSA_transcript": MSA_transcript_list,
                                             "projection": Projection_list,
                                             "Classification": Classification_list})
    return MSA_classifications

def download_gzipped_file(file_URL, name = "codon_alignments.fa.gz"):
    """Download gzipped file from URL into the current working directory.
    args:
        file_URL: URL where gzipped file can be found. URL should already contain
        the file name.
        name: The filename to give the downloaded gzipped file.
    """
    # disable verification warning
    disable_warnings(InsecureRequestWarning)
    # Request file
    web_response = requests.get(file_URL, verify=False)
    # download fasta file
    open(name, 'wb').write(web_response.content)
    # print message
    print(f'{name} downloaded from {file_URL}\n')

    return None

def select_sequence(fasta, fasta_keys, identifier, select = "QUERY"):
    """Function to retrieve sequence from identifier.
    args:
        fasta: fasta file, read from pyfastx.
        fasta_keys: the keys (i.e. identifiers) of the sequences in the fasta file.
        This way, it does not retrieve the keys again every time the function is called.
        idendifier: sequence identifier, does not need to be full identifier name. 
        'ENST00000667069' can also select the identifier 'ENST00000667069.MSH3.45'
        select: select the 'QUERY' or 'REFERENCE' sequence. Choose between these two.
    """
    # raise errors if arguments not correct
    if select not in ["QUERY", "REFERENCE"]:
        raise Exception("select argument must be 'QUERY' or 'REFERENCE'")
    if not type(select) is str:
        raise TypeError("select must be a string named 'QUERY' or 'REFERENCE'")
    if str(type(fasta)) != "<class 'Fasta'>":
        raise TypeError("fasta should be of of class Fasta")
    # filter ids for identifier
    ids_sel = fasta_keys.filter(fasta_keys % identifier)
    # return selected sequence
    for id in list(ids_sel):
        if select in id:
            return fasta[id]

def write_seq_to_fasta(seq, file, seq_header, append=False):
    """Function to write one sequence (retrieved by select_sequence() to a fasta file
    args:
        seq: the sequence to be read to a fasta file
        file: the filename, including file path, where the sequence should be stored.
        seq_header: the sequence header given to the sequence.
        append: if False, will write the sequence to a new fasta file (and delete already 
        existing files with the same name). If True, will append the sequence to an
        already existing fasta file.
    """
    # the fasta.seq method seems to not be working, so use the fasta.raw method
    raw_seq = str(seq.raw)
    # remove header, as a new sequence header wil be used
    raw_seq = raw_seq.replace(f'>{seq.name}\n', "")
    # remove spaces between the codons
    raw_seq = raw_seq.replace(" ", "")
    # write fasta file
    if append == False:
        mode = "w"
    elif append == True:
        mode = "r+"
    else:
        raise TypeError("append argument should be a boolean")
    with open(file, mode=mode, encoding="utf-8") as f:
        f.write(f'{seq_header}\n{raw_seq}')

def similarity_score(seq1, seq2):
    """Function to calculate the similarity score between two sequences,
    that were returned by the select_sequence() function in this module.

    args:
        The seq1 and seq2 sequences are two different sequences returned from
        the select_sequence() function.
    """
    # Extract raw string.
    raw_seq1, raw_seq2 = str(seq1.raw), str(seq2.raw)
    # remove header from raw string.
    raw_seq1, raw_seq2 = raw_seq1.replace(f'>{seq1.name}\n', ""), raw_seq2.replace(f'>{seq2.name}\n', "")
    # remove spaces between the codons
    raw_seq1, raw_seq2 = raw_seq1.replace(" ", ""), raw_seq2.replace(" ", "")

    # create aligner
    aligner = Align.PairwiseAligner()
    alignments = aligner.align(raw_seq1, raw_seq2)
    # pick the best alignment
    best_alignment = alignments[0]
    # calculate similarity score
    similarity = best_alignment.score / max(len(raw_seq1), len(raw_seq2))
    return similarity

def  QUERY_REFERENCE_similarity_scores(fa, projections):
    """Function to takes a list (or pd.series) of projections, and generates a dictionary of
    similarity scores. For every projection in the list, the similarity score between the QUERY and
    REFERENCE is calculated.

    args:
        fa: object of class fasta, read via the pyfastx.Fasta method using the pyfastx module.
        The fasta file should refer to a codonAlignments.fa.gz from 
        https://genome.senckenberg.de/download/TOGA/human_hg38_reference/.
        projections: list (or pd.series) of projections for which the similarity scores should be calculated.
    """
    similarity_scores = {}
    for projection in projections:
        
        seq_QUERY = select_sequence(fa, fa.keys(), identifier=projection, select="QUERY")
        seq_REFERENCE = select_sequence(fa, fa.keys(), identifier=projection, select="REFERENCE")

        sim_score = similarity_score(seq_QUERY, seq_REFERENCE)
        similarity_scores[projection] = sim_score
    
    return similarity_scores