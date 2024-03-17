#!/bin/bash

export subjects=(M9_4 M9_7 M10_1 M10_6 M11_4)
export thresholds=(5 4)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_intan_longterm_multiregion.py \
            --subject "$subject" \
            --sortdate 240317 \
            --threshold "$threshold"
    done
done
