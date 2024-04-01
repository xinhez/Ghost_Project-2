#!/bin/bash

export subjects=(L11_9 L16_8 L17_7)
export shanks=(-1 0 1 2 3 4)

for subject in "${subjects[@]}"
do  
    for shank in "${shanks[@]}"
    do
        python sort_intan.py \
            --subject "$subject" \
            --shank "$shank"
    done
done


export subjects=(L14_5)
export shanks=(0 1 2 3 4 5 6 7 -1)

for subject in "${subjects[@]}"
do 
    for shank in "${shanks[@]}"
    do
        python sort_intan.py \
            --subject "$subject" \
            --shank "$shank"
    done
done
