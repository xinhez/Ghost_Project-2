#!/bin/bash

export shanks=(-1) # (0 1 2 3 4)

export subjects=(M9_7 M10_1 M15_2 M15_5 M15_7 M16_1 M16_2 M17_5)
export thresholds=(4.5)
for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15 --sorted_region "region1"
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15 --sorted_region "region2"
        done
    done
done

export subjects=(M9_8 M10_6 M15_3 M16_6 M16_7 M17_2)
export thresholds=(5.5)
for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15 --sorted_region "region1"
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15 --sorted_region "region2"
        done
    done
done