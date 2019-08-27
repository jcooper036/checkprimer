#! /usr/bin/env python3
"""
Input: A folder that contains multiple csvs
Output: One csv that are all of them together
"""

import glob
import sys
import datetime


class Csv(object):

    def __init__(self, file, header=True):
        self.file = file
        self.header = header
        self.read_csv(self.header)
    
    def read_csv(self, header):
        """
        Input: csv file, header bool (optional)
        Gives the object the .data property
        """
        self.data = []
        with open(self.file, 'r') as f:
            for line in f:
                line = line.rstrip()
                if header:
                    self.header = line.split(',')
                if not header:
                    line = line.split(',')
                    self.data.append(line)
                header = False

def write_csv(csvlist, file, header=False):
    with open(file, 'w') as f:
        if header:
            f.write(','.join(header))
            f.write('\n')
        for ele in csvlist:
            f.write(','.join(ele))
            f.write('\n')

def combine_csvs(csvs):
    """
    """
    gotHeader = False
    master_list = []

    for csv in csvs:
        if not gotHeader:
            master_list.append(csvs[csv].header)
            gotHeader = True
        for ele in csvs[csv].data:
            master_list.append(ele)
    return master_list

# ID the folder
folder = sys.argv[1]

# use glob to find everything ending in .csv
csvs = {}
for csvFile in glob.glob(folder + "*/*_autoprimer_output.csv"):
    csvs[csvFile] = Csv(csvFile)

# get the contents of all those csvs together
output = combine_csvs(csvs)

# write teh output
now = datetime.datetime.now()
file = '/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer/autoprimer_combine_' + now.strftime("%Y-%m-%d.%H-%M") + '.csv'
write_csv(output, file)
