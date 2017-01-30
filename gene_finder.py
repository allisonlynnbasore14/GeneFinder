# -*- coding: utf-8 -*-
"""
YOUR HEADER COMMENT HERE

@author: Allison Basore

"""

import random
from amino_acids import aa, codons, aa_table   # you may find these useful
from load import load_seq


def shuffle_string(s):
    """Shuffles the characters in the input string
        NOTE: this is a helper function, you do not
        have to modify this in any way """
    return ''.join(random.sample(s, len(s)))

# YOU WILL START YOUR IMPLEMENTATION FROM HERE DOWN ###


def get_complement(nucleotide):

    """ Returns the complementary nucleotide

        nucleotide: a nucleotide (A, C, G, or T) represented as a string
        returns: the complementary nucleotide
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    """
    # TODO: implement this
    n = nucleotide
    if n == 'A':
        return 'T'
    if n =='G':
        return 'C'
    if n == 'C':
        return 'G'
    if n == 'T':
        return 'A'


def get_reverse_complement(dna):
    """ Computes the reverse complementary sequence of DNA for the specfied DNA
        sequence

        dna: a DNA sequence represented as a string
        returns: the reverse complementary DNA sequence represented as a string
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    """
    # TODO: implement this
    outputlist = ''
    rev_dna = dna[::-1]
    for i in rev_dna:
        n = get_complement(i)
        outputlist += n
    return outputlist

def rest_of_ORF(dna):
    """ Takes a DNA sequence that is assumed to begin with a start
        codon and returns the sequence up to but not including the
        first in frame stop codon.  If there is no in frame stop codon,
        returns the whole string.

        dna: a DNA sequence
        returns: the open reading frame represented as a string
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    """
    # TODO: implement this
    ORF = ''
    for i in range(0,len(dna),3):
        p = i +3
        bit = dna[i:p]
        if bit != 'TAG' and 'TAA' and 'TGA':
            ORF = ORF + dna[i:p]
        else:
            break
    return ORF



def find_all_ORFs_oneframe(dna):
    """ Finds all non-nested open reading frames in the given DNA
        sequence and returns them as a list.  This function should
        only find ORFs that are in the default frame of the sequence
        (i.e. they start on indices that are multiples of 3).
        By non-nested we mean that if an ORF occurs entirely within
        another ORF, it should not be included in the returned list of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    """
    # TODO: implement this
    #I am not 100% sure tha tthis function works properly. I am not sure if it will work with nested codons, if error, check
    output = []
    key = 'ATG' #setting the key for what I want to find
    for place in range(0,len(dna),3): #making slices of 3 place is start index
        bit_end = place+3 #bit_end is the stop index and changes with every new 'place' value
        bit = dna[place:bit_end] #takes a section of dna
        sub = False
        if bit == key:  #see if it matches the start codon
            newORF = rest_of_ORF(dna[place:len(dna)]) #if it is a start codon, it sends it to the rest of ORF function
            if len(output) > 0:
                for i in output:
                    if newORF in i:
                        sub = True
            if sub == False:
                output.append(newORF)
    return output


def find_all_ORFs(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence in
        all 3 possible frames and returns them as a list.  By non-nested we
        mean that if an ORF occurs entirely within another ORF and they are
        both in the same frame, it should not be included in the returned list
        of ORFs.

        dna: a DNA sequence
        returns: a list of non-nested ORFs

    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    """
    # TODO: implement this #here is the plan, start three frames and send them to the function: find_all_ORFs_oneframe
    frame1 = find_all_ORFs_oneframe(dna[0:len(dna)])
    frame2 = find_all_ORFs_oneframe(dna[1:len(dna)])
    frame3 = find_all_ORFs_oneframe(dna[2:len(dna)])
    output = []
    total = frame1
    for i in range(len(total)):
        sub = False
        string1 = total[i]
        for p in range(len(total)):
            string2 = total[p]
            if string1 in string2 and string1 != string2:
                sub = True
        if sub == False:
            output.append(string1)
    total2 = frame2
    for i in range(len(total2)):
        sub = False
        string1 = total2[i]
        for p in range(len(total)):
            string2 = total2[p]
            if string1 in string2 and string1 != string2:
                sub = True
        if sub == False:
            output.append(string1)
    total3 = frame3
    for i in range(len(total3)):
        sub = False
        string1 = total3[i]
        for p in range(len(total3)):
            string2 = total3[p]
            if string1 in string2 and string1 != string2:
                sub = True
        if sub == False:
            output.append(string1)
    return output


def find_all_ORFs_both_strands(dna):
    """ Finds all non-nested open reading frames in the given DNA sequence on both
        strands.

        dna: a DNA sequence
        returns: a list of non-nested ORFs
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    """
    # TODO: implement this
    #The plan is to call the reverse reverse compliment function and the find all ORF's fuction and make of list
    output = []
    strand1 = find_all_ORFs(dna)
    reverse = get_reverse_complement(dna)
    strand2 = find_all_ORFs(reverse)
    total = strand1 + strand2
    for i in range(len(total)):
        sub = False
        string1 = total[i]
        for p in range(len(total)):
            string2 = total[p]
            if string1 in string2 and string1 != string2:
                sub = True
        if sub == False:
            output.append(string1)
    return output


def longest_ORF(dna):
    """ Finds the longest ORF on both strands of the specified DNA and returns it
        as a string
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    """
    # TODO: implement this
    strands = find_all_ORFs_both_strands(dna) #gives the strand finder dna and gets backs the strans
    strands.sort(key=len) #sorts the strands from least to greatest
    return strands[len(strands)-1] #returns the last (greatest) one




def longest_ORF_noncoding(dna, num_trials):
    """ Computes the maximum length of the longest ORF over num_trials shuffles
        of the specfied DNA sequence

        dna: a DNA sequence
        num_trials: the number of random shuffles
        returns: the maximum length longest ORF """
    # TODO: implement this
    #this should return the length of the longest strand
    bankofstrands = []
    listofdna = list(dna)
    for i in range(num_trials):
        random.shuffle(listofdna)
        n = longest_ORF(listofdna)
        bankofstrands.append(n)
    bankofstrands.sort(key=len)
    return len(bankofstrands[len(bankofstrands)-1])

#is this how you test the previosu function?
#fin = open('X73525.fa')
#print(fin)
#longest_ORF_noncoding(fin,3)


def coding_strand_to_AA(dna):
    """ Computes the Protein encoded by a sequence of DNA.  This function
        does not check for start and stop codons (it assumes that the input
        DNA sequence represents an protein coding region).

        dna: a DNA sequence represented as a string
        returns: a string containing the sequence of amino acids encoded by the
                 the input DNA fragment

        >>> coding_strand_to_AA("ATGCGA")
        'MR'
        >>> coding_strand_to_AA("ATGCCCGCTTT")
        'MPA'
    """
    # TODO: implement this
    # This is going to return to the main genefinder function a string of the aminoacids from the dna given it.
    # The plan is to take the dna given, divide it into codones, pass these coddons to the table to look them up
    #    then store it all in a string that will be returned to the genefinder function. 
    stored_AAs = ''
    for i in range(0,len(dna),3):
        p = i +3
        bit = dna[i:p]
        if len(bit) < 3:
            break
        aminoA = aa_table[bit]
        stored_AAs = stored_AAs + aminoA
    return stored_AAs


def gene_finder(dna):
    """ Returns the amino acid sequences that are likely coded by the specified dna

        dna: a DNA sequence
        returns: a list of all amino acid sequences coded by the sequence dna.
    """
    # TODO: implement this
    # This function should call both the AA function and longest ORF function, I believe I should use the output from 
    #the longest ORF function into the amino acid function
    threshold = longest_ORF_noncoding(dna, 1500):
    larger_than_threshold = []
    aminolist = []
    findinglong = find_all_ORFs_both_strands(dna)
    if len(findinglong) > threshold:
        larger_than_threshold.append(findinglong)
    for i in larger_than_threshold:
        aminoA = coding_strand_to_AA(i)
        aminolist.append(aminoA)
    return aminolist

    #THE ABOVE FUNCTION I NOT TESTED. To do: test 'longestORF_noncoding', then test genefinder.


if __name__ == "__main__":
    import doctest
    doctest.run_docstring_examples(coding_strand_to_AA, globals(), verbose =True)