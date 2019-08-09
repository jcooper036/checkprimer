#! /usr/bin/env python3
"""
Input: A folder that contains multiple csvs
Output: One csv that are all of them together
"""


def read_csv(infile):
    csvList = []
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip()
            csvList.append(line)
    return csvList

def append_csv(csvList, outfile):
    with open(outfile, 'w+') as f:
        for line in csvList:
            f.write(line + '\n')

# use glob to find everything ending in .csv
csvs = []

# one by one, write them to a new csv. Do this s.t. only one csv is in memory at any given point

# write a blank new line to the master csv
master_csv = ''
with open(master_csv, 'w') as f:
    f.write('')

# loop over the csvs
for infile in csvs:
    csvList = read_csv(infile)
    append_csv(csvList, master_csv)

