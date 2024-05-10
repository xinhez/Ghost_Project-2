#!/bin/bash

export subjects=(D14_6) # D12_6 D13_4 D13_8 )

for subject in "${subjects[@]}"
do
    python sort_deepbrain_by_shank_drug.py --subject "$subject" --threshold 5.5
done