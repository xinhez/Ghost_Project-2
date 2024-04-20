#!/bin/bash

export subjects=(D13_4 D14_6 D12_6 D13_8 D14_4)
export thresholds=(3 3.5 4 4.5 5)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_deepbrain_by_channel.py \
            --subject "$subject" \
            --sortdate 240319 \
            --threshold "$threshold"
    done
done
