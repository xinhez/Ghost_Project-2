#!/bin/bash

export subjects=(M9_7 M10_6)
export thresholds=(3 5 4)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_intan_longterm_multiregion.py \
            --subject "$subject" \
            --savedate 240417 \
            --threshold "$threshold"
    done
done
