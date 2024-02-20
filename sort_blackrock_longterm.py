import matplotlib
matplotlib.use('Agg')

import argparse
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np 
import os
import spikeinterface.core as sc
import spikeinterface.extractors as se
import spikeinterface.curation as scu
import spikeinterface.preprocessing as spre 
import spikeinterface.sorters as ss
import sys 

from probeinterface import generate_multi_columns_probe
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes

sys.path.append('src')

from src.brpylib import NsxFile


# Constants 
n_s_per_min = 60
channel_indices = np.array([ 
    [ 5,  3,  1,  7,  9, 11], 
    [17, 15, 13, 19, 21, 23], 
    [29, 27, 25, 28, 26, 24], 
    [18, 20, 22, 16, 14, 12], 
    [ 6,  8, 10,  4,  2,  0],
])
shank_locations = np.array([[0, 0], [150, 200], [300, 400], [450, 200], [600, 0]])
sorter_parameters = {
    'mountainsort4': {
        'detect_sign': -1,
        'adjacency_radius': -1, 
        'freq_min': 300, 
        'freq_max': 3000,
        'filter': True,
        'whiten': True,  
        'clip_size': 50,
        'detect_threshold': 4,
        'detect_interval': 3, # 0.3ms 
    },
    'mountainsort5': {

    }
}
ms_before, ms_after = 2, 2

def get_args():
    """ Get command line arguments. """
    parser = argparse.ArgumentParser(description='Run Parameters')
    parser.add_argument(
        '--subject',
        type=str,
        help = 'subject name to sort',
    )
    parser.add_argument(
        '--data',
        type=str,
        help = 'root directory to host all files',
    )
    parser.add_argument(
        '--folder',
        type=str,
        help = 'directory to save all result',
    )
    parser.add_argument(
        '--sorter',
        type=str,
        default='mountainsort4',
        help = 'sorting algorithm',
    )
    parser.add_argument(
        '--nsx',
        type=str,
        default='ns6',
        help = 'file format to sort for',
    )
    parser.add_argument(
        '--savedate',
        type=str,
        default=datetime.datetime.today().strftime('%Y%m%d'),
        help = 'directory to save all result',
    )
    parser.add_argument(
        '--min_duration',
        type=int,
        default=n_s_per_min*10,
        help = 'minimum duration (s) to consider as valid data',
    )
    args = parser.parse_args()
    return args


def read_recordings(paths, min_duration, n_channel):
    traces = []
    files = []
    for path in paths:
        nsxfile = NsxFile(path)
        raw_data = nsxfile.getdata()
        nsxfile.close()
        trace = np.hstack(raw_data['data'])[:n_channel]
        sampling_frequency = raw_data['samp_per_s'] # assume the sampling frequency is always the same for the same file format.
        duration = sum([header['data_time_s'] for header in raw_data['data_headers']])
        if duration < min_duration:
            continue 
        traces.append(trace)
        files.append({
            'path': path,
            'sampling_frequency': sampling_frequency,
            'duration': duration,
        })
    traces = np.hstack(traces)
    return traces, files


def check_and_get_sampling_frequency(files):
    sampling_frequencies = [file['sampling_frequency'] for file in files]
    assert len(set(sampling_frequencies)) == 1 
    return sampling_frequencies[0]


def create_probe(channel_indices, shank_locations, savepath=None):
    plt.rcParams.update({'font.size': 10})
    n_shank = len(channel_indices)
    n_channel = channel_indices.size
    
    plt.figure(figsize=(20, 10))
    ax = plt.gca()

    probes = []
    for shank_channel_indices, shank_location in zip(channel_indices, shank_locations):
        probe = generate_multi_columns_probe(
            num_columns=2, num_contact_per_column=3, 
            xpitch=50, ypitch=50,
            contact_shapes='circle', contact_shape_params={'radius': 10}
        )
        probe.move(shank_location)
        probe.set_device_channel_indices(shank_channel_indices)

        plot_probe(probe, with_device_index=True, ax=ax)
        probes.append(probe)
        
    multi_shank_probe = combine_probes(probes)
    multi_shank_probe.set_device_channel_indices(channel_indices.flatten())

    plt.xlim(-100, 700)
    plt.ylim(-150, 550)
    plt.title(f'Probe - {n_channel}ch - {n_shank}shanks')
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    plt.close()
    return multi_shank_probe

def plot_unit(savepath):
    plt.savefig(savepath, bbox_to_inches='tight')
    plt.close()

if __name__ == '__main__':
    args = get_args()
    folder_root = f'{args.folder}/{args.subject}/{args.savedate}'
    print(f'Saving results to {folder_root}...')

    recording_folder = f'{folder_root}/recording'
    if not os.path.isfile(f'{recording_folder}/binary.json'):
        recording_paths = sorted(glob.glob(f'{args.data}/{args.subject}/**/*.{args.nsx}'))
        traces, files = read_recordings(recording_paths, args.min_duration, channel_indices.size)   
        sampling_frequency = check_and_get_sampling_frequency(files)
        recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
        recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
        recording = spre.common_reference(recording, reference='global', operator='median')
        recording.save(folder=recording_folder)
    recording = sc.load_extractor(recording_folder)
    probe = create_probe(channel_indices, shank_locations, f'{recording_folder}/probe.png')
    recording.set_probe(probe, in_place=True)
    print('\t...Preprocessed...')

    sorting_folder = f'{folder_root}/sorting'
    if not os.path.isfile(f'{sorting_folder}{os.sep}sorter_output{os.sep}firings.npz'):
        ss.run_sorter(
            sorter_name=args.sorter,
            recording=recording,
            output_folder = sorting_folder,
            remove_existing_folder=True,
            with_output=False,
            **sorter_parameters[args.sorter],
        )
    sorting = se.NpzSortingExtractor(f'{sorting_folder}{os.sep}sorter_output{os.sep}firings.npz')
    # spikeinterface https://github.com/SpikeInterface/spikeinterface/pull/1378
    sorting = scu.remove_excess_spikes(sorting, recording)
    print('\t...Sorted...')

    waveform_folder = f'{folder_root}/waveform'
    if not os.path.isfile(f'{waveform_folder}{os.sep}templates_average.npy'):
        sc.extract_waveforms(
            recording, sorting, 
            folder=waveform_folder,
            ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
            return_scaled=False,
            overwrite=True,
            use_relative_path=True,
        )
    print('\t...Waveform extracted...')
    waveform_extractor = sc.load_waveforms(
        folder=waveform_folder, with_recording=True, sorting=sorting
    )
    extremum_channels = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg') 

    # raise Exception
    # units_folder = f'{folder_root}/units'
    # os.makedirs(units_folder, exist_ok=True)
    # for unit_id in sorting.unit_ids:
    #     unit_plot_file = f'{units_folder}/{unit_id}.png'
    #     if not os.path.isfile(unit_plot_file):
    #         plot_unit(unit_id, recording, sorting, waveform_extractor, extremum_channels, savepath=unit_plot_file)