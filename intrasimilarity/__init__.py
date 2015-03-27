from __future__ import division

from subprocess import call
from pandas import DataFrame, Series
import os
from Bio import SeqIO

def cut_up_fasta(infasta, outfasta = None, chunk_size=1000, full_chunks = True):
    """
    returns non overlapping fasta chunks of length chunk_size
    if outfasta has a file name it writes it into the file
    if full_chunks is True, returns only the chunks that are "chunk_size" length (e.g. removes the terminal chunks)     
    """
    
    with open(infasta) as ff:
        seqs = []
        for record in SeqIO.parse(ff, "fasta"):
            for i in xrange(0,len(record),chunk_size):
                seqs += [ record[i:(i+chunk_size)] ]

    if full_chunks:
        seqs = [s for s in seqs if len(s) == chunk_size]

    for i,s in enumerate(seqs) :
        s.id = "seq_" + str(i)
        
    if outfasta:
        with open(outfasta, "w") as handle:
            SeqIO.write(seqs,handle,"fasta")
                
    return seqs

def NIC_similiarity(query_assembly, subject_assembly, chunk_size = 1000, identity_threshold = 95, length_threshold = 0.95):
    """
    a function computing similarity score between two assembly based on NICs (Near identicat Contigs). It chops up the query into chunks and counts how many of these chunks are NICs to the subject
    """
    chopped_up_query = "tmp.fasta"
    nb_chunks = len(cut_up_fasta(query_assembly, outfasta = chopped_up_query))
    nics = find_NICs(chopped_up_query, subject_assembly, identity_threshold, length_threshold)
    os.remove(chopped_up_query)
    return len(nics.keys())/nb_chunks

def find_NICs(query, subject, identity_threshold = 95, length_threshold = 0.95):
    """
    query is the path to the query
    subject is the path to the query
    keep only matches with at least identity_threshold nucleotide identity( default 99 )
    returns a dictionary with as key the seq-ids of the contigs that are highly similar over almost all there length to a part of an other content, and as value the contigs the map too
    requires blast+

    returns NICs (near identical contigs) of query in subject
    """
    
    blast_outp = "temp.tsv"
    word_size = "28"
    blast_db_files = [subject + ".nhr", subject + ".nin",  subject + ".nsq"]
    blast_db_cmd = ["makeblastdb" ,"-in", subject, "-dbtype", "nucl", "-out", subject]
    blast_cmd = ["blastn" , "-out",  blast_outp,  "-word_size", word_size , "-db", subject, "-query",  query, "-outfmt", "6 qseqid sseqid qlen slen pident length"]
    

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
    

def find_selfsimilars(assembly, identity_threshold = 95, length_threshold = 0.95):
    """
    assembly is the path to a fasta file to analyse
    keep only matches with at least identity_threshold nucleotide identity( default 99 )
    returns a dictionary with as key the seq-ids of the contigs that are highly similar over almost all there length to a part of an other content, and as value the contigs the map too
    requires blast+
    """

    return find_NICs(assembly,assembly,identity_threshold, length_threshold)
