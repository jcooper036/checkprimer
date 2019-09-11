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
GENOME = '/Volumes/i_bio/Crispr_F0_Screens/checkprimer/genome/GRCh38_latest_genomic.fasta'
GENOME = '/Users/jacob.cooper/resources/genomes/GRCh38_latest_genomic.fasta'
PRIMER_PER_SIDE = 4
END_BUFFER = 700
INSIDE_BUFFER = 150

import subprocess
import sys
import copy
import os

#########################
##### Classes
#########################

class Primer(object):

    def __init__(self, side, pen, seq, start, length, tm, gc, anyTH, endTH, hairpin, stab):
        self.side = side
        self.pen = pen
        self.seq = seq
        self.start = start
        self.length = length
        self.tm = tm
        self.gc = gc
        self.anyTH = anyTH
        self.endTH = endTH
        self.hairpin = hairpin
        self.stab = stab
    
    def output(self):
        out = "{},{},{},{},{}".format(self.seq, self.side, self.start, self.tm, self.gc)
        return out

#########################
##### Functions
#########################

def bash(command):
    """Makes running shell commands less verbose"""
    return subprocess.check_output(command, shell=True, executable='/bin/bash', universal_newlines=True)
    # return subprocess.check_output(command, shell=True, executable='/bin/bash', universal_newlines=True).decode('utf-8').strip()

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
    
        # add in the final entry
        if head and entry:
            fasta[head] = entry

    return fasta

def reverse_complement(sequence):
    """Reverse complement a DNA string"""
    rc_seq = ''
    rcDict = {
        'A' : 'T',
        'T' : 'A',
        'C' : 'G',
        'G' : 'C'
    }
    for let in sequence:
        rc_seq += rcDict[let.upper()]
    rc_seq = rc_seq[::-1]
    return rc_seq

def find_match(seq1, seq2):
    """
    Input: sequence 1 (primary sequence) and sequence 2 (the search)
    Returns: The start and end position of seq2 in seq1
    """

    found = False
    len2 = len(seq2)

    # try the sequence forward
    for idx, site in enumerate(seq1):
        if len(seq1) >= idx+len2:
            if seq1[idx:idx+len2] == seq2:
                found = True
                break
    if found:
        return (idx, idx+len2)
    
    # try the reverse complement
    seq2 = reverse_complement(seq2)
    for idx, site in enumerate(seq1):
        if len(seq1) >= idx+len2:
            if seq1[idx:idx+len2] == seq2:
                found = True
                break
    if found:
        return (idx, idx+len2)

def blast_primer(seq):
    """
    Input: primer sequence
    Returns True if primer is unique, false if the primer is not.
    """
    global GENOME
    if not os.path.exists(GENOME + '.nhr'):
        print("UPDATE: making BLAST database for " + GENOME)
        command = 'makeblastdb -in ' + GENOME + ' -parse_seqids -dbtype nucl'
        bash(command)
    command = 'blastn -task blastn-short -outfmt "6" -query <(echo ' + seq + ') -db ' + GENOME + ' | sort -k12 -n -r | head -5'
    out = str(bash(command)).split('\n')
    # see if there are hits greater than bit score of 36.2 (its the last entry in outfmt 6)
    examine = []
    for output in out:
        output = output.split('\t')
        if output[-1]:
            if float(output[-1]) > 36.2:
                examine.append(output)
    if len(examine) == 1:
        return True
    else:
        return False


def parse_primers(primer3out, side):
    primers = {}
    primer3out = primer3out.split('\n')
    global PRIMER_PER_SIDE
    count = 0
    max_blast = 8
    blastCount = 0
    tried_sequences = []
    reset = True
    print('Looking for {} primers'.format(side))
    for line in primer3out:
        if reset:
            pen = None
            seq = None
            start = None
            length = None
            gc = None
            tm = None
            anyTH = None
            endTH = None
            hairpin = None
            stab = None
            blastCount = 0
            reset = False

        if "PENALTY=" in line:
            pen = line.split('=')[1]
        if "SEQUENCE=" in line:
            seq = line.split('=')[1]
        if "," in line and 'EXPLAIN' not in line:
            # print(line)
            start = int(line.split('=')[1].split(',')[0])
            length = line.split('=')[1].split(',')[1]
        if "GC_PERCENT" in line:
            gc = line.split('=')[1]
        if "TM" in line:
            tm = line.split('=')[1]
        if "SELF_ANY" in line:
            anyTH = line.split('=')[1]
        if "SELF_END" in line:
            endTH = line.split('=')[1]
        if "HAIRPIN" in line:
            hairpin = line.split('=')[1]
        if "END_STABILITY" in line:
            stab = line.split('=')[1]

        # only save if the primer is more than 10bp away from all others in starting position
        if all([pen, seq, start, length, gc, tm, anyTH, endTH, hairpin, stab]):
            add = True
            
            if primers:
                overlapBuffer = 25
                starts = [int(primers[pr].start) for pr in primers]
                for s in starts:
                    if not (start >= s+overlapBuffer):
                        if not (start <= s-overlapBuffer):
                            # print('crowd catch')
                            add = False
            
            # make sure we haven't tried it before
            if add:
                if seq in tried_sequences:
                    add = False

            if (float(anyTH)+float(endTH)) > 8:
                add = False
            if add:
                if blastCount < max_blast:
                    print('Blasting {}, try {} / {}'.format(seq, blastCount+1, max_blast))
                    if not blast_primer(seq):
                        add = False
                    blastCount += 1
                else:
                    return primers
            if add:
                if side == 'left':
                    FR = 'F'
                if side == 'right':
                    FR = 'R'
                primers[FR+str(count+1)] = Primer(side, pen, seq, start, length, tm, gc, anyTH, endTH, hairpin, stab)
                count += 1
                reset = True
            # add it to the list of tried sequences REGUARDLESS
            tried_sequences.append(seq)

        if count == PRIMER_PER_SIDE:
            break
    return primers

def pick_primers(template, side, lowTm=60, highTm=65):
    """
    Input: Template sequence, side is left or right
    Returns: dictionary of primer objects, somewhere between 0 and 4
    """

    # settings variables
    sequence_id = 'blank'
    sequence_template = template
    # print(template) #@
    task = 'generic'
    if side == 'left':
        pick_left_primer = '1'
        pick_right_primer = '0'
    elif side == 'right':
        pick_left_primer = '0'
        pick_right_primer = '1'        
    
    pick_internal_oligo = '0'
    opt_size = '20'
    min_size = '19'
    max_size = '25'
    product_size_range = '75-150'
    explain_flag = '1'
    min_gc = '35'
    max_gc = '75'
    min_tm = '59'
    opt_tm = '63'
    max_tm = '67'
    primer_return_num = '200'

    # + 'PRIMER_PRODUCT_SIZE_RANGE=' + product_size_range + '\n'

    # string concat for settings
    primer3_settings = 'SEQUENCE_ID=' + sequence_id + '\n' + 'SEQUENCE_TEMPLATE=' + sequence_template + '\n' + 'PRIMER_TASK=' + task + '\n' + 'PRIMER_PICK_LEFT_PRIMER=' + pick_left_primer + '\n' + 'PRIMER_PICK_INTERNAL_OLIGO=' + pick_internal_oligo + '\n' + 'PRIMER_PICK_RIGHT_PRIMER=' + pick_right_primer + '\n' + 'PRIMER_OPT_SIZE=' + opt_size + '\n' + 'PRIMER_MIN_SIZE=' + min_size + '\n' + 'PRIMER_MAX_SIZE=' + max_size + '\n'  + 'PRIMER_EXPLAIN_FLAG=' + explain_flag + '\n' + 'PRIMER_INTERNAL_MIN_TM=' + min_tm + '\n' + 'PRIMER_INTERNAL_OPT_TM=' + opt_tm + '\n' + 'PRIMER_INTERNAL_MAX_TM=' + max_tm + '\n' + 'PRIMER_INTER_MIN_GC=' + min_gc + '\n' + 'PRIMER_INTERNAL_MAX_GC=' + max_gc + '\n' + 'PRIMER_NUM_RETURN=' + primer_return_num + '\n='

    # test for primer3
    return parse_primers(run_primer3(primer3_settings), side)

def output_primers(crisprs, gene, file):
    """
    Input: CRISPR dictionary (with the primers dictionaries inside it)
    Writes to a csv in the current directory
    """
    with open(file, 'w') as f:
        f.write('GENE,CRISPR,PRIMER,SEQUENCE,SIDE,START,TM,GC%\n')
        for cr in crisprs:
            for pr in crisprs[cr]['primers']:
                outstring = crisprs[cr]['primers'][pr].output()
                primerName = crisprs[cr]['primerTag'] + '_' + pr
                f.write("{},{},{},{}\n".format(crisprs[cr]['gene'], cr, primerName, outstring))


#########################
##### Run commands
#########################

# gather the sample and crispr targets
infile = sys.argv[1]
gene = infile.split('/')[-1].split('.fasta')[0]
fasta = read_fasta(infile)
crisprs = {}
for entry in fasta:
    if 'crispr' in entry:
        crseq = fasta[entry]
        entry = entry.split(' ')
        for info in entry:
            if 'crispr' in info:
                name = info.split(':')[1]
            if 'segment' in info:
                exon = info.split(':')[1]
        crisprs[name] = {
            'seq' : crseq,
            'name' : name,
            'segment' : exon,
            'gene' : gene,
            'primers' : {},
        }
        crisprs[name]['primerTag'] = '{}_{}_{}'.format(crisprs[name]['gene'], crisprs[name]['segment'], crisprs[name]['name'])

    else:
        seq = fasta[entry]

# terminal output for clarity
print('\n\n######################\nFinding primers for the CRISPR sites in {}'.format(gene))

# loop over each CRISPR molecule
for cr in crisprs:

    # terminal output for clarity
    print('\nLooking for primers for {}'.format(cr))

    # find the start and stop positions
    crisprs[cr]['start'], crisprs[cr]['stop'] = find_match(seq, crisprs[cr]['seq'])

    # find all the possible primers

    # left primers
    # first if there is pleantly of room
    end_buffer = copy.deepcopy(END_BUFFER)
    inside_buffer = copy.deepcopy(INSIDE_BUFFER)
    if (crisprs[cr]['start'] - end_buffer) > 0:
        
        # template is -end_buffer - -inside_buffer from crispr start site
        template = seq[crisprs[cr]['start']-end_buffer:crisprs[cr]['start']-inside_buffer]
        temp = pick_primers(template, 'left')
        for primer in temp:
            temp[primer].start = crisprs[cr]['start'] - inside_buffer - temp[primer].start
            crisprs[cr]['primers'][primer] = copy.deepcopy(temp[primer])
    
    ## make sure one of the left primers is close enough for sequencing
    check = False
    for primer in crisprs[cr]['primers']:
        p_start = crisprs[cr]['primers'][primer].start
        c_start = crisprs[cr]['start']
        if abs(p_start-c_start) < 500:
            check = True
    if not check:    
        end_buffer = 499
        
    # right primers
    # first if there is pleantly of room
    if (crisprs[cr]['stop'] + end_buffer) < len(seq):
        
         # template is +inside_buffer to +end_buffer from crispr start site
        template = seq[crisprs[cr]['stop']+inside_buffer:crisprs[cr]['stop']+end_buffer]
        temp = pick_primers(template, 'right')
        for primer in temp:
            temp[primer].start = crisprs[cr]['stop'] + inside_buffer + temp[primer].start
            crisprs[cr]['primers'][primer] = copy.deepcopy(temp[primer])

    # write to the terminal
    print(crisprs[cr]['name'])
    for primer in crisprs[cr]['primers']:
        name = primer
        primer_seq = crisprs[cr]['primers'][primer].seq
        start = crisprs[cr]['primers'][primer].start
        print(name, primer_seq, start)

# output the seqeunces to a csv in the current directory
file = infile.split('.fasta')[0] + "_autoprimer_output.csv"
output_primers(crisprs, gene, file)