# Checkprimers
This program will use primer3 to design check primers around a given sequence of DNA. In this specific case, it is being used to find check primers around a sgRNA site.

## Criteria
- primers need to be non-overlapping.
- they should amplify a region 500-1000bp
- they should be at minimum 50bp from the check region, and at maximum 400bp
- should all be a consistent Tm

## General Steps
1) Read an input file that contains:
    - Maybe a single sequence and target to design around
    - Maybe many sequences?

2) For each input, run a primer3 command
    - By specifing that we only pick either left or right primers at a time, we shoud be able to narrow the region quite easily.

3) If we don't find a primer, iterate until we do or we run out of sequence to search.
    - Question here of how much parameters should / can be expanded for the search.

4) Output the primers. This should probably be in a csv.
    - Each primer should just go with the gene name
    - Gets an ID#, F or R
    - ID number could correspond to gene, or just be gene name. Ask the users.