#!/bin/bash

export subjects=(1_5 5_7 6_2 6_7 7_2 8_6) 

for subject in "${subjects[@]}"
do  
    python sort_blackrock_longterm.py \
        --subject "$subject" \
        --data data/raw/LongTerm \
        --folder data/processed/LongTerm \
        --nsx ns4
done
