import glob 
import os 
import shutil 

date = '240319'
output_folder = f'data/processed/{date}-dp-resort'

thresholds = {
    'D12_6': 4.0,
    'D13_4': 4.0,
    'D13_8': 4.0,
    'D14_6': 4.0,
}

for subject, threshold in thresholds.items():
    for src in (glob.glob(f'data/processed/{subject}/{date}/traces-curation-units/shank*-{threshold}')):
        subject =  src.split('/')[2]
        param = src.split('/')[-1]
        tgt = f'{output_folder}/{subject}/{param}'
        print(src)
        print(tgt)
        shutil.copytree(src, tgt)
        print()