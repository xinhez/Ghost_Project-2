#!/bin/bash

export subjects=(D13_4 D14_6)
export thresholds=(4 5)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_deepbrain.py \
            --subject "$subject" \
            --sortdate 240319 \
            --threshold "$threshold"
    done
done
