import glob
import os 
import shutil 

sortdate = '240310'


output_folder = f'data/{sortdate}'
os.makedirs(output_folder, exist_ok=True)
for path in sorted(glob.glob(f'data/processed/MultiRegion/**/*{sortdate}*/**/units*')):
    subject = path.split('/')[3]
    session = path.split('/')[4]
    region = path.split('/')[-2]
    sorting_param = path.split('/')[-1]

    tgt = f'{output_folder}/{subject}/{region}/{sorting_param}'
    print(path)
    print(tgt)
    print()