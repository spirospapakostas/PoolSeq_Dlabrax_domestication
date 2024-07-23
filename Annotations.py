#!/usr/bin/env python
# coding: utf-8

# In[13]:


# Download the BioMart annotations
# install gff3_parser from within jupyter 
# gff3_parser is a simple python package to parse gff3 (Generic Feature Format) files into pandas dataframes.
# more info at: https://github.com/McClain-Thiel/gff3_parser


# In[11]:


import sys
get_ipython().system('{sys.executable} -m pip install gff3-parser')
import gff3_parser

filepath = "Dicentrarchus_labrax.seabass_V1.0.105.gff3"
dicentrachus_biomart_annotations_v1 = gff3_parser.parse_gff3(filepath, verbose = True, parse_attributes = True)


# In[ ]:


# import warnings
warnings.filterwarnings("ignore")  # disable warnings

LG_ID = "HG916841.1" # HG916841.1 is an exmaple LG_ID
SNP_position= 19222880 # an eaxmple position of SNP


dicentrachus_biomart_annotations= dicentrachus_biomart_annotations_v1[(dicentrachus_biomart_annotations_v1["Seqid"]==LG_ID)]
dicentrachus_biomart_annotations["Start"] = dicentrachus_biomart_annotations["Start"].astype('int')
dicentrachus_biomart_annotations["End"] = dicentrachus_biomart_annotations["End"].astype('int')
a= dicentrachus_biomart_annotations[(dicentrachus_biomart_annotations["Start"]>  SNP_position  -100000) &
                                    (dicentrachus_biomart_annotations["End"]<  SNP_position  +100000)]
a

