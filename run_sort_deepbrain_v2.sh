#!/bin/bash

export subjects=(D13_8 D12_6 D13_4 D14_6)

for subject in "${subjects[@]}"
do
    python sort_deepbrain_v2.py --subject "$subject"
done