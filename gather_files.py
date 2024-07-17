import glob
import os 
import shutil


threshold = 5.5
sorted_duration = 15

# output_folder = f'data/processed/units{threshold}-by-shank'
# os.makedirs(output_folder, exist_ok=True)

# for plot_file in (glob.glob(f'data/processed/**/shank*/region*/units{threshold}-{sorted_duration}min/*.png')):
#     subject = plot_file.split('/')[2]
#     shank = plot_file.split('/')[3]
#     region = plot_file.split('/')[4]
#     unit_id = plot_file.split('/')[-1].split('.')[0]
#     tgt = f'{output_folder}/{subject}/{region}/{shank}-unit{unit_id}.png'

#     os.makedirs(f'{output_folder}/{subject}/{region}/', exist_ok=True)
#     shutil.copy2(plot_file, tgt)
# print(len((glob.glob(f'data/processed/**/shank*/region*/units{threshold}-{sorted_duration}min/*.png'))))



output_folder = f'data/processed/units{threshold}'
os.makedirs(output_folder, exist_ok=True)

for plot_file in (glob.glob(f'data/processed/**/all/region*/units{threshold}-{sorted_duration}min/*.png')):
    subject = plot_file.split('/')[2]
    region = plot_file.split('/')[4]
    unit_id = plot_file.split('/')[-1].split('.')[0]
    tgt = f'{output_folder}/{subject}/{region}/unit{unit_id}.png'

    os.makedirs(f'{output_folder}/{subject}/{region}/', exist_ok=True)
    shutil.copy2(plot_file, tgt)
print(len((glob.glob(f'data/processed/**/all/region*/units{threshold}-{sorted_duration}min/*.png'))))