#!/bin/bash

export subjects=(L14_5)
export shanks=(2)
export thresholds=(5.0)

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

export subjects=(L11_9)
export shanks=(0 1 2 4)
export thresholds=(5.0)

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

export subjects=(L17_7)
export shanks=(0 1)
export thresholds=(5.0)

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