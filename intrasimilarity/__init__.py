from __future__ import division

from subprocess import call
from pandas import DataFrame, Series
import os


def find_selfsimilars(assembly, identity_threshold = 99, length_threshold = 0.90):
    """
    assembly is the path to a fasta file to analyse
    keep only matches with at least identity_threshold nucleotide identity( default 99 )
    returns a dictionary with as key the seq-ids of the contigs that are highly similar over almost all there length to a part of an other content, and as value the contigs the map too
    requires blast+
    """
    
    blast_outp = "temp.tsv"
    word_size = "28"
    blast_db_files = [assembly + ".nhr", assembly + ".nin",  assembly + ".nsq"]
    blast_db_cmd = ["makeblastdb" ,"-in", assembly, "-dbtype", "nucl", "-out", assembly]
    blast_cmd = ["blastn" , "-out",  blast_outp,  "-word_size", word_size , "-db", assembly, "-query",  assembly, "-outfmt", "6 qseqid sseqid qlen slen pident length"]
    

    blastdb_return = call(blast_db_cmd)
    blast_return = call(blast_cmd)

    # open blast data
    blast_data = DataFrame.from_csv(blast_outp, sep = "\t", header=None, index_col=None)
    blast_data.columns = Series(["qseqid", "sseqid", "qlen", "slen", "pident", "length"])

    #remove hit on the same contig
    blast_data = blast_data.loc[blast_data['qseqid'] != blast_data['sseqid']]

    #filter out "low" similarity
    blast_data = blast_data.loc[blast_data['pident'] > identity_threshold]

    #keep only matches where the length of the query is a significant proportion of the match (according to length_threshold)
    blast_data = blast_data.loc[blast_data['length']/blast_data['qlen'] > length_threshold]
    
    os.remove(blast_outp)
    for f in blast_db_files:
        os.remove(f)

    #returns a dictionary with as key the seq-ids of the contigs that are highly similar over almost all there length to a part of an other content, and as value the contigs the map too
    return {q : list(set(blast_data.loc[blast_data['qseqid'] == q]['sseqid'])) for q in blast_data['qseqid']}
