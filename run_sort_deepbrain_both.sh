#!/bin/bash

export subjects=(D12_6 D13_4 D13_8 D14_6)
export thresholds=(5.5 5.0 4.5 6.0)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do
        python sort_deepbrain_both.py --subject "$subject" --threshold "$threshold"
    done
done