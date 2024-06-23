#!/bin/bash

export shanks=(-1)
export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
export thresholds=(4.5)
for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --do_sorting 0 --plot_traces 0
        done
    done
done

# export shanks=(-1)
# export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
# export thresholds=(4.5)
# for threshold in "${thresholds[@]}"
# do
#     for subject in "${subjects[@]}"
#     do  
#         for shank in "${shanks[@]}"
#         do
#             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --do_sorting 0 --plot_traces 1
#         done
#     done
# done
########################################################################################

export shanks=(0 1 2 3 4)
export thresholds=(4.5)
export regionOneSubjects=(M16_1 M15_7 M15_5 M15_3 M15_2 M10_1)
for threshold in "${thresholds[@]}"
do
    for subject in "${regionOneSubjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region1" --sorted_duration 15
        done
    done
done

export regionTwoSubjects=(M9_7 M9_8 M10_6 M15_3 M15_5 M15_7 M16_2)
for threshold in "${thresholds[@]}"
do
    for subject in "${regionTwoSubjects[@]}"
    do
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region2" --sorted_duration 15
        done
    done
done

export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15
        done
    done
done

########################################################################################

export shanks=(0 1 2 3 4)
export thresholds=(3.5)
export regionOneSubjects=(M16_1 M15_7 M15_5 M15_3 M15_2 M10_1)
for threshold in "${thresholds[@]}"
do
    for subject in "${regionOneSubjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region1" --sorted_duration 15
        done
    done
done

export regionTwoSubjects=(M9_7 M9_8 M10_6 M15_3 M15_5 M15_7 M16_2)
for threshold in "${thresholds[@]}"
do
    for subject in "${regionTwoSubjects[@]}"
    do
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region2" --sorted_duration 15
        done
    done
done

export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15
        done
    done
done

########################################################################################

export shanks=(0 1 2 3 4)
export thresholds=(5.5)
export regionOneSubjects=(M16_1 M15_7 M15_5 M15_3 M15_2 M10_1)
for threshold in "${thresholds[@]}"
do
    for subject in "${regionOneSubjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region1" --sorted_duration 15
        done
    done
done

export regionTwoSubjects=(M9_7 M9_8 M10_6 M15_3 M15_5 M15_7 M16_2)
for threshold in "${thresholds[@]}"
do
    for subject in "${regionTwoSubjects[@]}"
    do
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region2" --sorted_duration 15
        done
    done
done

export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
for threshold in "${thresholds[@]}"
do
    for subject in "${subjects[@]}"
    do  
        for shank in "${shanks[@]}"
        do
            python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_duration 15
        done
    done
done

########################################################################################


# export shanks=(-1)

# # export subjects=(M9_7 M9_8 M10_1 M10_6 M15_2 M15_3 M15_5 M15_7 M16_1 M16_2 M16_6 M16_7 M17_2 M17_5)
# # export thresholds=(4.5)
# # for threshold in "${thresholds[@]}"
# # do
# #     for subject in "${subjects[@]}"
# #     do  
# #         for shank in "${shanks[@]}"
# #         do
# #             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --do_sorting 0
# #         done
# #     done
# # done

# export regionOneSubjects=(M16_1 M15_7 M15_5 M15_3 M15_2 M10_1)
# export thresholds=(4.5)
# for threshold in "${thresholds[@]}"
# do
#     for subject in "${regionOneSubjects[@]}"
#     do  
#         for shank in "${shanks[@]}"
#         do
#             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region1" --sorted_duration 15
#         done
#     done
# done

# export regionTwoSubjects=(M9_7 M9_8 M10_6 M15_3 M15_5 M15_7 M16_2)
# export thresholds=(4.5)
# for threshold in "${thresholds[@]}"
# do
#     for subject in "${regionTwoSubjects[@]}"
#     do
#         for shank in "${shanks[@]}"
#         do
#             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region2" --sorted_duration 15
#         done
#     done
# done


# export regionOneSubjects=(M16_1 M15_7 M15_5 M15_3 M15_2 M10_1)
# export thresholds=(5.5)
# for threshold in "${thresholds[@]}"
# do
#     for subject in "${regionOneSubjects[@]}"
#     do  
#         for shank in "${shanks[@]}"
#         do
#             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region1" --sorted_duration 15
#         done
#     done
# done

# export regionTwoSubjects=(M9_7 M9_8 M10_6 M15_3 M15_5 M15_7 M16_2)
# export thresholds=(5.5)
# for threshold in "${thresholds[@]}"
# do
#     for subject in "${regionTwoSubjects[@]}"
#     do
#         for shank in "${shanks[@]}"
#         do
#             python sort_intan.py --subject "$subject" --shank "$shank" --threshold "$threshold" --sorted_region "region2" --sorted_duration 15
#         done
#     done
# done

# # export doubleRegionSubjects=(M15_3 M15_5 M15_7)
# export thresholds=(4.5)