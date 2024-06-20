#!/bin/bash

export subjects=(L14_5)
export shanks=(0 1 2 3 4 5 6 7 -1)
export thresholds=(3)

for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do 
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold"
        done
    done
done

export subjects=(L16_8 L17_7 L11_9) 
export shanks=(0 1 2 3 4 -1)
export thresholds=(3)

for threshold in "${thresholds[@]}"
do
   for subject in "${subjects[@]}"
   do  
       for shank in "${shanks[@]}"
       do
           python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold"
   	done
   done
done