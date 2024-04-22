import glob 
import os 
import shutil 

date = '240319'
output_folder = f'data/processed/{date}-dp'

for src in (glob.glob(f'data/processed/**/{date}/units*')):
    subject =  src.split('/')[2]
    param = src.split('/')[-1]
    tgt = f'{output_folder}/{subject}/{param}'
    print(src)
    print(tgt)
    shutil.copytree(src, tgt)
    print()