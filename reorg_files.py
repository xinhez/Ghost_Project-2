import os 
import shutil

root = 'data/raw/Lateral_tosort'
outroot = 'data/raw'
dates = sorted(os.listdir(root))

for date in dates:
    sessions = sorted(os.listdir(f'{root}/{date}'))
    for session in sessions:
        mouse = '_'.join(session.split('_')[:2])
        print('date', date)
        print('mouse', mouse)
        src = f'{root}/{date}/{session}'
        if os.path.isdir(src):
            tgt = f'{outroot}/{mouse}/{date}/{session}'
            print('src', repr(src))
            print('tgt', repr(tgt))
            print()
            shutil.move(src, tgt)
    print()