import os 
import shutil 

thresholds = [3.5, 4.5, 5.5, 6.5]
mice = ['M15_2', 'M15_3', 'M15_5', 'M15_7', 'M16_1']
regions = ['CA1', 'M1']

output_folder = 'sorted'

for threshold in thresholds:
    for mouse in mice:
        for region in regions:
            src = f'data/processed/{mouse}/{region}/units-{threshold}'
            tgt = f'{output_folder}/{mouse}/{region}/units-{threshold}'
            print(src)
            print(tgt)
            print()
            shutil.copytree(src, tgt)