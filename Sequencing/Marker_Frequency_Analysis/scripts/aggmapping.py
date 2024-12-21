import csv
import copy
import ast
import re
import os
import numpy as np
import pandas as pd
import seaborn as sns
import dask.dataframe as dd

from Bio import SeqIO
from Bio.Seq import Seq

from matplotlib import pyplot as plt

def cigar_length(cigar_string):
    """
    Calculate the total length of the sequence represented by a CIGAR string.
    
    Parameters:
    cigar_string (str): The CIGAR string to be quantified.
    
    Returns:
    int: The total length of the sequence.
    """
    import re

    # Regular expression to match CIGAR operations
    cigar_re = re.compile(r'(\d+)([MIDNSHP=X])')
    
    total_length = 0

    # Iterate over all matches of the CIGAR operations
    for length, operation in cigar_re.findall(cigar_string):
        length = int(length)
        if operation in 'MD':
            # M (match/mismatch), D (deletion) contribute to length
            total_length += length

    return total_length

def readpositionsfromsam(samfilepath):
    read_position_df = []
    with open(samfilepath, "r") as samfile:
        for line in samfile:
            if line[0] == "@":
                next(samfile)
            else:
                splitline = line.split("\t")
                position = splitline[3]
                qual = int(splitline[4])
                cigar = splitline[5]
                cigar_len = cigar_length(cigar)
                mid_position = int(position) + (cigar_len//2)
                read_position_df.append([mid_position,qual])                
    read_position_df = pd.DataFrame(read_position_df,columns=["Ref Position","MAPping Quality"])
    return read_position_df

sampath = snakemake.input["alignmentpath"]
samquantpath = snakemake.output["samquantpath"]
read_position_df_final = []
for chunk_file in os.listdir(sampath):
    chunkpath = sampath + "/" + chunk_file
    read_position_df = readpositionsfromsam(chunkpath)
    read_position_df_final.append(read_position_df)
read_position_df_final = pd.concat(read_position_df_final)
read_position_df_final.to_pickle(samquantpath)



    
