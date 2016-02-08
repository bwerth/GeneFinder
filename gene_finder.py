# -*- coding: utf-8 -*-
"""
This is a gene finding program that determines regions of the Salmonella bacterium's DNA that code for proteins

@author: Bryan Werth

"""
import math
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
    """
    Returns the complement of a nucleotide given an input of A, C, G, B.
    Anything else and the function returns None
    """
    """
    >>> get_complement('A')
    'T'
    >>> get_complement('C')
    'G'
    >>> get_complement('G')
    'C'
    >>> get_complement('B')
    None
    """
    if nucleotide == 'A':
        return 'T'
    elif nucleotide == 'T':
        return 'A'
    elif nucleotide == 'C':
        return 'G'
    elif nucleotide == 'G':
        return 'C'
    else:
        print 'This is not a nucleotide'
        return None

def get_reverse_complement(dna):
    """
    This function returns the reverse complementary DNA sequence for the input DNA sequence. 
    """
    """
    >>> get_reverse_complement("ATGCCCGCTTT")
    'AAAGCGGGCAT'
    >>> get_reverse_complement("CCGCGTTCA")
    'TGAACGCGG'
    >>> get_reverse_complement("")
    ''
    """
    n = len(dna)
    reverse_complement_dna = ''
    while n>0:
        complement = get_complement(dna[n-1])
        reverse_complement_dna = reverse_complement_dna + complement
        n = n - 1
    return reverse_complement_dna

def rest_of_ORF(dna):
    """
    Takes an input sequence of DNA that is assumed to begin with a start codon, 
    and returns the snippet of DNA from the beginning of the string up to, but 
    not including, the first in frame stop codon.  If there is no in frame stop 
    codon, the whole string is returned.
    """
    """
    >>> rest_of_ORF("ATGTGAA")
    'ATG'
    >>> rest_of_ORF("ATGAGATAGG")
    'ATGAGA'
    >>> rest_of_ORF("TAG")
    ''
    """
    i=0
    while i<len(dna):
        codon = dna[i:i+3]
        if codon == 'TAG' or codon == 'TAA' or codon == 'TGA':
            return dna[:i]
        else:
            i = i + 3
    return dna

def find_all_ORFs_oneframe(dna):
    """
    Finds all non-nested open reading frames in the given DNA
    sequence and returns them as a list.  This function should
    only find ORFs that are in the default frame of the sequence
    (i.e. they start on indices that are multiples of 3).
    By non-nested we mean that if an ORF occurs entirely within
    another ORF, it should not be included in the returned list of ORFs.
    """
    """            
    >>> find_all_ORFs_oneframe("ATGCATGAATGTAGATAGATGTGCCC")
    ['ATGCATGAATGTAGA', 'ATGTGCCC']
    >>> find_all_ORFs_oneframe("CAAATGAATGTAGATAGATGTGCCC")
    ['ATGAATGTAG','ATGTGCCC']
    >>> find_all_ORFs_oneframe("")
    []
    >>> find_all_ORFs_oneframe("AATGTAG")
    []
    """
    i=0
    all_ORFs = []
    while i<len(dna):
        codon = dna[i:i+3]
        if codon == 'ATG':
            all_ORFs.append(rest_of_ORF(dna[i:]))
            i=i+len(rest_of_ORF(dna[i:]))
        else:
            i = i + 3
    return all_ORFs


def find_all_ORFs(dna):
    """
    Finds all non-nested open reading frames in the given DNA sequence in
    all 3 possible frames and returns them as a list.  By non-nested we
    mean that if an ORF occurs entirely within another ORF and they are
    both in the same frame, it should not be included in the returned list
    of ORFs.
    """
    """
    >>> find_all_ORFs("ATGCATGAATGTAG")
    ['ATGCATGAATGTAG', 'ATGAATGTAG', 'ATG']
    >>> find_all_ORFs("AATGTAG")
    ['ATG']
    >>> find_all_ORFs("")
    []
    >>> find_all_ORFs("AAAAA")
    []
    """
    String0 = find_all_ORFs_oneframe(dna)
    String1 = find_all_ORFs_oneframe(dna[1:])
    String2 = find_all_ORFs_oneframe(dna[2:])
    return String0 + String1 + String2

def find_all_ORFs_both_strands(dna):
    """
    This finds ORFs on both the original DNA sequence and its reverse complement.
    """
    """
    >>> find_all_ORFs_both_strands("ATGCGAATGTAGCATCAAA")
    ['ATGCGAATG', 'ATGCTACATTCGCAT']
    >>> find_all_ORFs_both_strands("")
    []
    >>> find_all_ORFs_both_strands("GATGTA")
    ["ATG"]
    >>> find_all_ORFs_both_strands("AAAAA")
    []
    """
    Strand1 = find_all_ORFs(dna)
    Strand2_RC = get_reverse_complement(dna)
    Strand2 = find_all_ORFs(Strand2_RC)
    return Strand1 + Strand2

def longest_ORF(dna):
    """
    Finds the longest open reading frame on either strand of the DNA
    """
    """
    >>> longest_ORF("ATGCGAATGTAGCATCAAA")
    'ATGCTACATTCGCAT'
    >>> longest_ORF("ATGTAGGATAAAGTA")
    'ATGAAA'
    >>> longest_ORF("")
    ''
    >>> longest_ORF("AAA")
    ''
    """
    ORF_List = find_all_ORFs_both_strands(dna)
    max_index = len(ORF_List)
    n = 0
    long_ORF = ''
    while n<max_index:
        if len(long_ORF)<len(ORF_List[n]):
            long_ORF = ORF_List[n]
        n = n + 1
    return long_ORF

def random_shuffle(dna):
    """
    Randomly shuffles a string
    """
    l = list(dna)
    random.shuffle(l)
    return ''.join(l)

def longest_ORF_noncoding(dna,num_of_trials):
    """
    Takes as input a dna strand and an integer number of trials to be executed. It randomly shuffles the 
    dna strand and computes the longest ORF overall over the course of the given number of trials
    """
    n = -1
    testing_ORF = ''
    long_ORF = ''
    while n<num_of_trials:
        dna = random_shuffle(dna)
        if len(longest_ORF(dna))>len(long_ORF):
            long_ORF = longest_ORF(dna)
        n = n + 1
    return len(long_ORF)

def coding_strand_to_AA(dna):
    """
    This function converts from a string containing a DNA sequence to a sequence of amino acids.
    """
    """
    >>> coding_strand_to_AA("ATGCGA")
    'MR'
    >>> coding_strand_to_AA("")
    ''
    """
    i = 0
    AA = ''
    while i<len(dna):
        codon = dna[i:i+3]
        if len(codon)<3:
            return AA
        else:
            aa = aa_table[codon]
            AA = AA + aa
            i = i + 3
    return AA

def gene_finder(dna):
    """
    Gene Finder takes a dna strand as input and returns a list of amino acid sequences that are
    longer than the threshold, which is taken using the longest_ORF_noncoding function.
    """
    i = 0
    threshold = longest_ORF_noncoding(dna, 1500)
    print threshold
    ORF_List = find_all_ORFs_both_strands(dna)
    AA_List = []
    while i<len(ORF_List):
        if len(ORF_List[i]) > threshold:
            AA_List.append(coding_strand_to_AA(ORF_List[i]))
        i = i + 1
    return AA_List

if __name__ == "__main__":
    import doctest
    doctest.testmod()

from load import load_seq
dna = load_seq("./data/X73525.fa")
print coding_strand_to_AA('ATGAAATAG')
