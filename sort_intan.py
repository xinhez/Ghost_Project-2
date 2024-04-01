import matplotlib 
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000

import argparse
import glob
import sys

sys.path.append('src')
from src.facts import probe_designs

sorter_parameters = {
    'detect_sign': -1,
    'adjacency_radius': -1, 
    'freq_min': 300, 
    'freq_max': 3000,
    'filter': False,
    'whiten': True,  
    'clip_size': 50,
    'num_workers': 8,
}

def get_args():
    parser = argparse.ArgumentParser(description='Run Parameters')
    parser.add_argument(
        '--subject',
        type=str,
        help='subject name to sort',
    )
    parser.add_argument(
        '--shank', 
        type=int,
        default=-1,
        help='shank number to sort, -1 for all shanks'
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=3,
        help='sorting detect threshold',
    )
    parser.add_argument(
        '--min_duration',
        type=int,
        default=10,
        help='minimum duration (min) to consider as valid data',
    )
    args = parser.parse_args()
    return args


def main(args):
    output_root = f'data/processed/{args.subject}/{"all" if args.shank < 0 else f"shank{args.shank}"}'
    segment_paths = sorted(glob.glob(f'data/raw/{args.subject}/**/**'))
    print(f'Saveing to {output_root}')
    print(f'Sorting {args.subject} with {len(segment_paths)} segment(s):')
    for segment_index, segment_path in enumerate(segment_paths):
        print(f'    [{segment_index+1:3d}] {segment_path}')
    
    sorter_parameters['detect_threshold'] = args.threshold

    if 'multiregion' in probe_designs[args.subject]:
        from sort_intan_multiregion import sort 
        sort(args, output_root, segment_paths, sorter_parameters)
    elif 'singleregion' in  probe_designs[args.subject]:
        from sort_intan_singleregion import sort 
        sort(args, output_root, segment_paths, sorter_parameters)
    else:
        raise Exception("Not implemented!!!")

if __name__ == '__main__':
    args = get_args()
    main(args)
