#! /usr/bin/env python3

import sys

def translate_sequence(rna_sequence, genetic_code):
    rna_seqs = rna_sequence.upper()
    protein = ""
    if len(rna_seqs) >= 3:
        """for i in range (0, len(rna_seqs), 3):
            codon = rna_seqs[i:i + 3]
            protein+= genetic_code[codon]"""
        if len(rna_seqs)%3 == 1:
            rna_seqs = rna_seqs[:-1]
        for i in range (0, len(rna_seqs), 3):
            codon = rna_seqs[i:i + 3]
            protein+= genetic_code[codon]
        if '*' in str(protein):
            sep = '*'
            protein = protein.split(sep, 1)[0]
            protein =  protein.replace('*', '')
    return protein
    pass


def get_all_translations(rna_sequence, genetic_code):
    rna_seqs = rna_sequence.upper()
    aa_list = []
    start_pos = 0

    def translate(start_pos,rna_seqs,genetic_code):
        proteins = ''
        for i in range(start_pos, len(rna_seqs), 3):
            codon = rna_seqs[i:i + 3]
            if codon in ['UAG', 'UAA', 'UGA'] or len(codon) != 3:
                break
            else: proteins += genetic_code[codon]
        return proteins

    while start_pos < len(rna_seqs):
        start_codon = rna_seqs[start_pos:start_pos + 3]
        if start_codon == 'AUG':
            translation = translate(start_pos, rna_seqs, genetic_code)
            aa_list.append(translation)
        start_pos += 1
    return aa_list
#    if rna_seqs.count('AUG') == 0:
#        protein = []
#    if rna_seqs.count('AUG') == 1:
#        protein = ['M']
#    if rna_seqs.count('AUG') == 2:
#        rna_seqs = rna_seqs.split('AUG')
#        if rna_seqs.count('') == 1:
#            del rna_seqs[0]
#            for section in rna_seqs:
#                codon = rna_seqs[section:section + 3]
#                protein+= genetic_code[codon]


        #if rna_seqs.count('') == 0:
        #    rna_seqs
        #rna_seqs = ''.join(rna_seqs)
        #for i in range (0, len(rna_seqs), 3):
        #    codon = rna_seqs[i:i + 3]
        #    protein+= genetic_code[codon]
        #    protein = ''.join(protein)
    #rna_seqs = rna_seqs.split(sep, 1)[0]
    #rna_seqs_list = rna_seqs.split('AUG')
    #sep = 'AUG'
    #rna_seqs_list = [e+sep for e in rna_seqs_list(sep) if e]
    #if len(rna_seqs_list) >= 3:
    """for i in range (0, len(rna_seqs), 3):
            codon = rna_seqs[i:i + 3]
            protein+= genetic_code[codon]"""
        #if len(rna_seqs)%3 == 1:
            #rna_seqs = rna_seqs[:-1]
    #    for name in rna_seqs_list:
         # (0, len(rna_seqs_list), 3):
    #        codon = rna_seqs_list[name:name + 3]
    #        protein+= genetic_code[codon]
    #        protein = ['M'+name for name in protein(name) if e]
    #    for i in range(0,len(protein)):
    #        protein[i] = str(protein[i])
    #return protein
    pass
 #if '*' in str(protein):
            #sep = '*'
            #protein = protein.split(sep, 1)[0]
            #protein =  protein.replace('*', '')
    #sep = 'M'
    #protein = protein.split(sep, 1)[0] 
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """

def get_reverse(sequence):
    seqs = list(sequence.upper())
    rev_seqs = seqs[::-1]
    reversed = ("".join(rev_seqs))
    return reversed
    pass #get reverse

def get_complement(sequence):
    sequence = list(sequence.upper())
    complement = ""
    for i in sequence:
        if i == "A":
            complement += "U"
        if i == "U":
            complement += "A"
        if i == "C":
            complement += "G"
        if i == "G":
            complement += "C"
    return complement
    pass #get complement

def reverse_and_complement(sequence):
    seqs = list(sequence.upper())
    rev_seqs = seqs[::-1]
    reversed = ("".join(rev_seqs))
    complement = ""
    for i in reversed:
        if i == "A":
            complement += "U"
        if i == "U":
            complement += "A"
        if i == "C":
            complement += "G"
        if i == "G":
            complement += "C"
    return complement
    pass #reverse and complement

def get_longest_peptide(rna_sequence, genetic_code):
    rna_seqs = rna_sequence.upper()
    aa_list = []
    start_pos = 0

    def translate(start_pos,rna_seqs,genetic_code):
        proteins = ''
        for i in range(start_pos, len(rna_seqs), 3):
            codon = rna_seqs[i:i + 3]
            if codon in ['UAG', 'UAA', 'UGA'] or len(codon) != 3:
                break
            else: proteins += genetic_code[codon]
        return proteins

    def reverse_and_complement(rna_seqs):
        rna_seqs = list(rna_seqs.upper())
        rev_seqs = rna_seqs[::-1]
        reversed = ("".join(rev_seqs))
        complement = ""
        for i in reversed:
            if i == "A":
                complement += "U"
            if i == "U":
                complement += "A"
            if i == "C":
                complement += "G"
            if i == "G":
                complement += "C"
        return complement

    rev_and_comp = reverse_and_complement(rna_seqs)
#    translation = translate(start_pos, rna_seqs, genetic_code)
#   translation_back = translate(start_pos, rev_and_comp, genetic_code)

    while start_pos < len(rna_seqs):
        start_codon = rna_seqs[start_pos:start_pos + 3]
        if start_codon == 'AUG':
            translation = translate(start_pos, rna_seqs, genetic_code)
            aa_list.append(translation)
        start_codon = rev_and_comp[start_pos:start_pos + 3]
        if start_codon == 'AUG':
            translation_back = translate(start_pos, rev_and_comp, genetic_code)
            aa_list.append(translation_back)
        start_pos += 1

#    while start_pos < len(rev_and_comp):
#        start_codon = rev_and_comp[start_pos:start_pos + 3]
#        if start_codon == 'AUG':
#            translation_back = translate(start_pos, rev_and_comp, genetic_code)
#            aa_list.append(translation_back)
#        start_pos += 1

    check = 'M'
    aa_list_start =  [idx for idx in aa_list if idx[0].lower().startswith(check.lower())]
    longest_aa = ''
    if aa_list_start == []:
        longest_aa == ''
    if aa_list_start != []:
        longest_aa = max((aa_list_start), key=len)
    return longest_aa
    pass
"""Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """



if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")

