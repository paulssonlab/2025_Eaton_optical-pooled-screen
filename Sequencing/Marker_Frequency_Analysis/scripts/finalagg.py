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

sampaths = snakemake.input["samquantpaths"]
samquantpathagg = snakemake.output["samquantpathagg"]

all_sam_df = []
for sampath in sampaths:
    samfilename = sampath.split("/")[-1]
    fullsamplename = samfilename[:-4]
    sam_df = pd.read_pickle(sampath)
    sam_df["Full Sample Name"] = fullsamplename
    sam_df["Full Sample Name"] = sam_df["Full Sample Name"].astype("category")
    sam_df["Ref Position"] = sam_df["Ref Position"].astype("uint32")
    sam_df["MAPping Quality"] = sam_df["MAPping Quality"].astype("uint16")
    all_sam_df.append(sam_df)

all_categories = set.union(*[set(item["Full Sample Name"].cat.categories) for item in all_sam_df])
for i in range(len(all_sam_df)):
    all_sam_df[i]["Full Sample Name"] = all_sam_df[i]["Full Sample Name"].cat.set_categories(all_categories)
    
all_sam_df = pd.concat(all_sam_df).reset_index(drop=True)
all_sam_df.to_pickle(samquantpathagg)