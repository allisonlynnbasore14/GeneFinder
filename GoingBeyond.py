# -*- coding: utf-8 -*-
"""
Created on 1/29/2017

@author: Allison Basore
"""

from os import path
#First loading the dna sequence that is known to code for Nitrogenase
from load import load_nitrogenase_seq
#This is how you get the meta-genome
from load import load_metagenome

#Goal of this program is to determine which of these DNA tuples match the nitrogen sequence
#Find the longest matching one, dosent have to match exactly

def loading_dna():
    metagenome = load_metagenome()#So, this actually is a list of tubles with the name of a seqquence and then the sequence.
    pass

def substring_checkc():
    """Returns the parts of the string that match"""
    nitrogenase = load_nitrogenase_seq()
    metagenome = load_metagenome()#So, this actually is a list of tubles with the name of a seqquence and then the sequence.
    
    for i in metagenome:
        print(i)
    pass

def finding_semi_match():
    """calls other functions to find a near match to the nitrogen gene, returns the tuple of the dna that mathces"""

    pass