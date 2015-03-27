from __future__ import division

import subprocess
from pandas import DataFrame, Series
import os
from Bio import SeqIO
import numpy as np
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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
    blast_cmd = ["blastn" , "-out",  blast_outp,  "-word_size", word_size , "-db", subject, "-query",  query, "-perc_identity", str(identity_threshold), "-outfmt", "6 qseqid sseqid qlen slen pident length"]
    
    with open("/dev/null") as null:
        blastdb_return = subprocess.call(blast_db_cmd, stdout=null)
    blast_return = subprocess.call(blast_cmd)

    # open blast data
    blast_data = DataFrame.from_csv(blast_outp, sep = "\t", header=None, index_col=None)
    blast_data.columns = Series(["qseqid", "sseqid", "qlen", "slen", "pident", "length"])

    blast_data = blast_data.loc[blast_data['pident'] > identity_threshold]

    
    #remove hit on the same contig
    blast_data = blast_data.loc[blast_data['qseqid'] != blast_data['sseqid']]


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

def test_run_with_gage_assemblies():
    """
    a little function running the NIC similarity for the Staph Assemblies of GAGE
    one can plot the output of this function with plot_heatmap(test_run_with_gage_assemblies())
    """
    # test data to be found here : http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Assembly.tgz

    data_folder = "/home/murumbi/repos/AssemblyValidation/Assembly/"
    scaffs =[ "ABySS2" , "Allpaths-LG" , "MSR-CA" , "SOAPdenovo" , "ABySS" , "Bambus2" , "SGA" , "Velvet" ]
    asses = {k: data_folder + k + "/genome.ctg.fasta" for k in scaffs}
    
    return  compare_assemblies(asses)



def compare_assemblies(assemblies):
    """
    compares a set of assemblies:
    assemblies is a dictionary with names of the assemblies as keys and fasta-files of the assemblies as values
    """
    similarities = {}
    for scaff_name, scaff in assemblies.iteritems():
        similarities[scaff_name] = {} 
        for subject_name, subject in assemblies.iteritems():
            print scaff_name, "vs", subject_name
            similarities[scaff_name][subject_name] = NIC_similiarity(scaff, subject, chunk_size = 100)

    similars =  DataFrame.from_dict(similarities)
    return similars

def plot_heatmap(similars):
    # row_clustering
    pairwise_dists = distance.squareform(distance.pdist(similars))
    row_clusters = sch.linkage(pairwise_dists,method='complete')
    den_rows = sch.dendrogram(row_clusters,color_threshold=np.inf,no_plot=True)

    # col_clustering
    col_pairwise_dists = distance.squareform(distance.pdist(similars.T))
    col_clusters = sch.linkage(col_pairwise_dists,method='complete')
    den_cols = sch.dendrogram(col_clusters,color_threshold=np.inf,no_plot=True)
    
    # reorder dataframe
    similars = similars.ix[den_rows['leaves']]
    similars = similars[den_cols['leaves']]
    
    fig = plt.figure()
    heatmapGS = gridspec.GridSpec(2,2,wspace=0.0,hspace=0.0,width_ratios=[0.25,1],height_ratios=[0.25,1])
    col_denAX = fig.add_subplot(heatmapGS[0,1])
    col_denD = sch.dendrogram(col_clusters,color_threshold=np.inf)
#    clean_axis(col_denAX)

    ### row dendrogram ###
    row_denAX = fig.add_subplot(heatmapGS[1,0])
    row_denD = sch.dendrogram(row_clusters,color_threshold=np.inf,orientation='right')
#    clean_axis(row_denAX)
    
    heatmapAX = fig.add_subplot(heatmapGS[1,1])
    axi = heatmapAX.imshow(similars,interpolation='nearest',aspect='auto',origin='lower',cmap=plt.cm.RdBu)
#    clean_axis(heatmapAX)
    heatmapAX.set_xticks(np.arange(0,len(similars.columns)))
    heatmapAX.set_xticklabels(similars.columns, rotation='vertical')

    heatmapAX.set_yticks(np.arange(0,len(similars.index)))
    heatmapAX.yaxis.set_ticks_position('right')
    heatmapAX.set_yticklabels(similars.index)

    ### legend

    scale_cbGSSS = gridspec.GridSpecFromSubplotSpec(1,2,subplot_spec=heatmapGS[0,0],wspace=0.0,hspace=0.0)
    scale_cbAX = fig.add_subplot(scale_cbGSSS[0,0]) # colorbar for scale in upper left corner
    cb = fig.colorbar(axi,scale_cbAX) # note that we tell colorbar to use the scale_cbAX axis
    cb.set_label('NIC-similarity')
    cb.ax.yaxis.set_ticks_position('left') # move ticks to left side of colorbar to avoid problems with tight_layout
    cb.ax.yaxis.set_label_position('left') # move label to left side of colorbar to avoid problems with tight_layout
    cb.outline.set_linewidth(0)

    fig.tight_layout()
    plt.show()

def clean_axis(ax):
    """Remove ticks, tick labels, and frame from axis"""
    ax.get_xaxis().set_ticks([])
    ax.get_yaxis().set_ticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
