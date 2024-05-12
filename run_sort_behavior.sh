#!/bin/bash

export subjects=(M16_1 M15_7 M15_5 M15_3 M15_2)
export thresholds=(5.5 6.5 4.5 3.5)
export run_sorts=(False True)

for threshold in ${thresholds[@]}
do
    for run_sort in ${run_sorts[@]}
    do 
        for subject in "${subjects[@]}"
        do  
            python sort_behavior.py --subject "$subject" --threshold $threshold --run_sort $run_sort
        done
    done
done
