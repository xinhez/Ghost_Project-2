import os 
import shutil

root = 'data/raw/LongTerm_tosort'
outroot = 'data/raw/LongTerm'
dates = sorted(os.listdir(root))

for date in dates:
    mice = sorted(os.listdir(f'{root}/{date}'))
    print(date)
    print(mice)
    for mouse in mice:
        src = f'{root}/{date}/{mouse}'
        if os.path.isdir(src):
            tgt = f'{outroot}/{mouse}/{date}'
            print('src', src)
            print('tgt', tgt)
            print()
            shutil.move(src, tgt)
    print()