#!/bin/bash

export thresholds=(3 3.5 4 4.5 5)

for threshold in "${thresholds[@]}"
do
    python sort_deepbrain_by_shank.py --threshold "$threshold"
done