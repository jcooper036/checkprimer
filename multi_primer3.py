#! /usr/bin/env python3

"""
The goal of this sript is to take an input
    - file of sequence and their CRISPR primers
And return PCR check primers for each one. The check primers need to:
    - Produce products that are 500bp - 1000bp
    - ideally all have the same (or very similar) Tms
    - all sequences need 2 forward and 2 reverse primers
    - non-overlapping
    - at least 100bp between the primer and the CRISPR target site

"""

"""
Here is what the primer3 example input looks like
SEQUENCE_ID=example
SEQUENCE_TEMPLATE=tattggtgaagcctcaggtagtgcagaatatgaaacttcaggatccagtgggcatgctactggtagtgctgccggccttacaggcattatggtggcaaagtcgacagagttta
PRIMER_TASK=generic
PRIMER_PICK_LEFT_PRIMER=1
PRIMER_PICK_INTERNAL_OLIGO=0
PRIMER_PICK_RIGHT_PRIMER=1
PRIMER_OPT_SIZE=20
PRIMER_MIN_SIZE=18
PRIMER_MAX_SIZE=22
PRIMER_PRODUCT_SIZE_RANGE=75-150
PRIMER_EXPLAIN_FLAG=1
=
"""

"""
Ideally we organize the information so it can be streamed into primer3
"""

import subprocess

def bash(command):
    """Makes running shell commands less verbose"""
    return subprocess.check_output(command, shell=True, universal_newlines=True)

def run_primer3(primer3_settings):
    """Appends the necessary quotes to the primer3 settings, runs primer3, returns the shell output"""
    primer3_settings = '"' + primer3_settings + '"'
    out = bash('printf ' + primer3_settings + ' | ~/primer3/src/primer3_core')
    return out

def read_fasta(file):
    """
    Input: fasta formatted file
    Returns: dictionary of fasta file
    """
    
    fasta = {}
    with open(file, 'r') as f:
        # set some flags and variables
        header = False
        entry = ''
        head = ''
        
        # loop over the lines
        for line in f:
            line = line.rstrip()

            ## check if a new header
            if '>' in line:

                # add the recorded entry
                if head and entry:
                    fasta[head] = entry
                # set the header to true
                header = True

            # reset the header
            if header:
                entry = ''
                head = line.split('>')[1]
                header = False
            # otherwize keep adding to the entry
            else:
                entry += line
    return fasta

def find_match(seq1, seq2):
    """
    Input: sequence 1 (primary sequence) and sequence 2 (the search)
    Returns: The start and end position of seq2 in seq1
    """
    
    len2 = len(seq2)
    for idx, site in enumerate(seq1):
        if len(seq1) >= idx+len2:
            if seq1[idx:idx+len2] == seq2:
                return (idx, idx+len2)


def pick_primers(template, side, lowTm=60, highTm=65):
    """
    Input: Template sequence, side is left or right
    Returns: dictionary of primer objects, somewhere between 0 and 4
    """

    # settings variables
    sequence_id = 'blank'
    sequence_template = template
    task = 'generic'
    if side == 'left':
        pick_left_primer = '1'
        pick_right_primer = '0'
    elif side == 'right':
        pick_left_primer = '0'
        pick_right_primer = '1'        
    
    pick_internal_oligo = '0'
    opt_size = '20'
    min_size = '18'
    max_size = '22'
    product_size_range = '75-150'
    explain_flag = '1'

    # string concat for settings
    primer3_settings = 'SEQUENCE_ID=' + sequence_id + '\n' + 'SEQUENCE_TEMPLATE=' + sequence_template + '\n' + 'PRIMER_TASK=' + task + '\n' + 'PRIMER_PICK_LEFT_PRIMER=' + pick_left_primer + '\n' + 'PRIMER_PICK_INTERNAL_OLIGO=' + pick_internal_oligo + '\n' + 'PRIMER_PICK_RIGHT_PRIMER=' + pick_right_primer + '\n' + 'PRIMER_OPT_SIZE=' + opt_size + '\n' + 'PRIMER_MIN_SIZE=' + min_size + '\n' + 'PRIMER_MAX_SIZE=' + max_size + '\n' + 'PRIMER_PRODUCT_SIZE_RANGE=' + product_size_range + '\n' + 'PRIMER_EXPLAIN_FLAG=' + explain_flag + '\n='

    # test for primer3
    print(run_primer3(primer3_settings))


# going to need a regex statement to find the primers, since they are PRIMER_SIDE_X, where X is a number. 


# gather the sample and crispr targets
fasta = read_fasta('multi_primer3/examples/SCGN/SCGN.fasta')
crisprs = {}
for entry in fasta:
    if len(fasta[entry]) < 25 and len(fasta[entry]) > 20:
        crisprs[entry] = {
            'seq' : fasta[entry],
            'name' : entry
        }
    else:
        seq = fasta[entry]

# find the start and end position of each CRISPR in the sequences
for cr in crisprs:
    crisprs[cr]['start'], crisprs[cr]['stop'] = find_match(seq, crisprs[cr]['seq'])

# pick left primers
for cr in crisprs:
    # first if there is pleantly of room
    if (crisprs[cr]['start'] - 400) > 0:
        template = seq[crisprs[cr]['start']-400:crisprs[cr]['start']-100]
        pick_primers(template, 'left')

# pick right primers
for cr in crisprs:
    # first if there is pleantly of room
    if (crisprs[cr]['stop'] + 400) < len(seq):
        template = seq[crisprs[cr]['stop']+100:crisprs[cr]['stop']+400]
        pick_primers(template, 'right')