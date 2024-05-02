#!/bin/bash

export subjects=(1_5 5_7 6_2 6_7 7_2) # 
export shanks=(0 1 2 3 4)
export thresholds=(3 4)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_blackrock.py \
                --subject "$subject" \
                --shank "$shank" \
		        --threshold "$threshold"
    	done
    done
done

export subjects=(5_7 6_2 6_7 7_2) # 1_5
export shanks=(-1)
export thresholds=(3 4)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_blackrock.py \
                --subject "$subject" \
                --shank "$shank" \
		        --threshold "$threshold"
    	done
    done
done