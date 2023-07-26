import pandas as pd
import numpy as np
import os
import json
import re
import requests

url_dict = {"uniHEART":"http://heart.ifr.fidt.top:61010/api/system/search/ugtd"}
meta_dict = {
    "uniHEART":['cell_ID','donor_ID','donor_gender','donor_age','original_name','organ','region','subregion','sample_status','seq_tech','cell_type','if_patient','donor_status','treatment','ethnicity','Ref','MCT','develop_stage']
}
ref_gene = pd.read_table("total_gene_list_43878.txt")
# default headers
headers={
    "Accept": "application/json, text/plain, */*",
    "Accept-Encoding": "gzip, deflate",
    "Cache-Control": "no-cache",
    "Content-Length": "255",
    "Content-Type": "application/json",
    "Host": "heart.ifr.fidt.top:61010",
    "Origin": "http://heart.ifr.fidt.top:61010",
    "Pragma": "no-cache",
    "Proxy-Connection": "keep-alive",
    "Referer": "http://heart.ifr.fidt.top:61010/",
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Safari/537.36"
}


def connect_ugtd(atlas):
    global url
    global meta_columns
    if atlas in ["uniHEART"]:
        url = url_dict[atlas]
        meta_columns = meta_dict[atlas]
    else:
        print('Current version of ECAUGT only availble for "uniHEART".')

def list_metadata_name():
    return meta_columns
        
def check_metadata_name(col_name):
    if col_name in meta_columns:
        print("The column " + col_name + " is in the database.")
    else:
        print("The column " + col_name + " is NOT in the database. Please more details of the available metadata from the paper.")
        
def check_gene_name(gene_name):
    if gene_name[0:4] == "ENSG":
        print("Please input a gene symbol instead of the ensembl ID.")
        if gene_name in ref_gene['Ensembl ID'].values:
            print("You may want to check the gene: " + ref_gene.loc[ref_gene['Ensembl ID']==gene_name,"Gene Symbol"].values[0])
    elif gene_name in ref_gene['Gene Symbol'].values:
        print("The gene " + gene_name + " is in the database.")
    else:
        print("The gene " + gene_name + " is NOT in the database. Please more details of the available genes from the paper.")
        
def _standardize_seq(seq):
    # delete spaces in the string
    seq = seq.strip(" ")
    seq = re.sub('\s*!\s*','!',seq)
    seq = re.sub('\s*\(\s*','(',seq)
    seq = re.sub('\s*\)\s*',')',seq)
    seq = re.sub('\s*&&\s*','&&',seq)
    seq = re.sub('\s*\|\|\s*','||',seq)   
    seq = re.sub('\s*==\s*','==',seq)
    seq = re.sub('\s*<>\s*','<>',seq)
    seq = re.sub('\s*>\s*','>',seq)
    seq = re.sub('\s*<\s*','<',seq)
    seq = re.sub('\s*>=\s*','>=',seq)
    seq = re.sub('\s*<=\s*','<=',seq)
    return seq

def get_data(conditions, columns_to_get=["cid"], do_transfer = False):
    # limit the number of the columns to get
    if len(columns_to_get)>2000:
        print("If you want to get more than 2000 columns, we suggest to use graphic interface to download the data.")
        return -1
    for column in columns_to_get:
        if (column in meta_columns) or (column in ref_gene['Gene Symbol'].values) or (column == 'cid'):
            pass
        else:
            print("Column {} not in the database.".format(column))
            return -1
    # columns_to_get should include 'cid'
    if not "cid" in columns_to_get:
        columns_to_get.extend(["cid"])
        
    # standardize condition string
    conditions_std = _standardize_seq(conditions)
    
    # construct POST data
    data={
        "conditions": conditions_std,
        "columns": columns_to_get
    }
    data = json.dumps(data)
    recieved  = requests.post(url=url, data=data, headers=headers).json()
    if recieved["ret"]:
        if do_transfer:
            df = pd.DataFrame(recieved['result'])
            return df.loc[:,columns_to_get]
        else:
            return recieved['result']
    else:
        print("Error: Get data failed.")
        return recieved["msg"]
    return recieved