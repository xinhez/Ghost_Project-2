#!/bin/bash

export subjects=(M16_1 M15_7 M15_5 M15_3 M15_2)
export thresholds=(5.5 4.5 3.5)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        python sort_behavior.py --subject "$subject" --threshold "$threshold"
    done
done
