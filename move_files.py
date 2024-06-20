import glob 
import os 
import shutil


output_folder = f'data/traces'
for path in glob.glob('data/processed/**/all/region*/traces'):
    # print(path)
    subject = path.split('/')[2]
    region = path.split('/')[4]

    output_path = f'{output_folder}/{subject}/{region}'
    print(path)
    print(output_path)
    print()
    shutil.copytree(path, output_path)