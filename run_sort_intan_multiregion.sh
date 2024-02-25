#!/bin/bash

export subjects=(M10_5 M10_8 M11_4 M9_4 M9_7 M10_1)
export thresholds=(4)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_intan_multiregion.py \
            --subject "$subject" \
            --sortdate 240225 \
            --threshold "$threshold"
    done
done