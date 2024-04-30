import glob
import os 
import shutil 

sortdate = '240417'

########################################################################################################################
# Single Session
########################################################################################################################
# output_folder = f'data/{sortdate}-multi'
# os.makedirs(output_folder, exist_ok=True)
# for path in sorted(glob.glob(f'data/processed/MultiRegion/**/*{sortdate}*/**/traces.png')):
#     subject = path.split('/')[3]
#     session = path.split('/')[4]
#     region = path.split('/')[-2]
#     sorting_param = path.split('/')[-1]

#     tgt = f'{output_folder}/{subject}/{region}.png'
#     os.makedirs(f'{output_folder}/{subject}', exist_ok=True)
#     if os.path.isfile(tgt):
#         tgt = tgt.replace('.png', '-2.png')
#     print(path)
#     print(tgt)
#     print()
#     shutil.copy2(path, tgt)


# output_folder = f'data/{sortdate}-multi'
# os.makedirs(output_folder, exist_ok=True)
# for path in sorted(glob.glob(f'data/processed/MultiRegion/**/*{sortdate}*/**/units*.0')):
#     subject = path.split('/')[3]
#     session = path.split('/')[4]
#     region = path.split('/')[-2]
#     sorting_param = path.split('/')[-1]

#     tgt = f'{output_folder}/{subject}/{region}/{sorting_param}'
#     os.makedirs(f'{output_folder}/{subject}', exist_ok=True)
#     print(path)
#     print(tgt)
#     print()
#     if os.path.isdir(tgt):
#         tgt = tgt + '-2'
#     shutil.copytree(path, tgt)



# output_folder = f'data/{sortdate}-multi'
# os.makedirs(output_folder, exist_ok=True)
# for path in sorted(glob.glob(f'data/processed/MultiRegion/**/*{sortdate}*/**/session_traces')):
#     subject = path.split('/')[3]
#     session = path.split('/')[4]
#     region = path.split('/')[-2]
#     sorting_param = path.split('/')[-1]

#     tgt = f'{output_folder}/{subject}/{region}/{sorting_param}'
#     os.makedirs(f'{output_folder}/{subject}', exist_ok=True)
#     print(path)
#     print(tgt)
#     print()
#     if os.path.isdir(tgt):
#         tgt = tgt + '-2'
#     shutil.copytree(path, tgt)

########################################################################################################################
# Single Session
########################################################################################################################



########################################################################################################################
# Multi Session
########################################################################################################################

output_folder = f'data/{sortdate}-multi-longterm'
os.makedirs(output_folder, exist_ok=True)
for path in sorted(glob.glob(f'data/processed/longterm/**/{sortdate}/**/units*')):
    subject = path.split('/')[3]
    session = path.split('/')[4]
    region = path.split('/')[-2]
    sorting_param = path.split('/')[-1]

    tgt = f'{output_folder}/{subject}/{region}/{sorting_param}'
    os.makedirs(f'{output_folder}/{subject}', exist_ok=True)
    print(path)
    print(tgt)
    print()
    if os.path.isdir(tgt):
        tgt = tgt + '-2'
    shutil.copytree(path, tgt)

########################################################################################################################
# Multi Session
########################################################################################################################