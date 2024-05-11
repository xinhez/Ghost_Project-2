import glob 
import os 
import shutil 

date = '240510'
threshold = 5.5
output_folder = f'data/processed/{date}-experiment-both-{threshold}'


# for subject in ['D12_6', 'D13_4', 'D13_8', 'D14_6']:
#     for tag, group in zip(['saline', 'drug'], ['240506', '240507']):
#         for src in (glob.glob(f'data/processed/{subject}/{group}/units-5.5')):
#             tgt = f'{output_folder}/{subject}/{tag}'
#             print(src)
#             print(tgt)
#             shutil.copytree(src, tgt)
#             print()

for subject in ['D12_6', 'D13_4', 'D13_8', 'D14_6']:
    for src in (glob.glob(f'data/processed/experiment/{subject}/units-{threshold}')):
        tgt = f'{output_folder}/{subject}/'
        print(src)
        print(tgt)
        shutil.copytree(src, tgt)
        print()