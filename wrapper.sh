#!/bin/bash

# wrapper for executing all multiprimer design in the /Volumes/i_bio/Crispr_F0_Screens/0-Genes\ for\ design/Genes_for_autoprimer directory

#/Volumes/i_bio/Crispr_F0_Screens/0-Genes\ for\ design/
DIR=/Volumes/i_bio/Crispr_F0_Screens/0-Genes_for_design/Genes_for_autoprimer/

for f in ${DIR}*/*.fasta
do
    if /Volumes/i_bio/Crispr_F0_Screens/checkprimer/multi_primer3.py ${f}; then
        echo 'pass'
    else
        echo 'FAILURE: '${f}
    fi
done    

/Volumes/i_bio/Crispr_F0_Screens/checkprimer/concat_csv.py ${DIR}0-Genes\ for\ design/Genes_for_autoprimer/