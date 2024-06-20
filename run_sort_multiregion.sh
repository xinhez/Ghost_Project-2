#!/bin/bash

export shanks=(-1)
export thresholds=(4)

export regionTwoSubjects=(M9_7 M9_8 M10_6 M15_3 M15_5 M15_7 M16_2)

for threshold in "${thresholds[@]}"
do
    for subject in "${regionTwoSubjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --region "region2"
        done
    done
done

export regionOneSubjects=(M10_1 M15_2 M15_3 M15_5 M15_7 M16_1)
export thresholds=(5)

for threshold in "${thresholds[@]}"
do
    for subject in "${regionOneSubjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --region "region1"
        done
    done
done

# export doubleRegionSubjects=(M15_3 M15_5 M15_7)
# export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_5)

# for threshold in "${thresholds[@]}"
# do
#     for subject in "${subjects[@]}"
#     do  
#         for shank in "${shanks[@]}"
#         do
#             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold"
#         done
#     done
# done
