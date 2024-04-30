#!/bin/bash

export subjects=(M9_7 M9_8 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
export thresholds=(4 5)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_intan_multiregion.py \
            --subject "$subject" \
            --sortdate 240428 \
            --threshold "$threshold"
    done
done
