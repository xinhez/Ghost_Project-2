import matplotlib 
matplotlib.use('Agg')
matplotlib.rcParams['agg.path.chunksize'] = 10000

import anndata as ad
import argparse
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
import os 
import pandas as pd
import re
import scanpy
import shutil
import spikeinterface.core as sc 
import spikeinterface.curation as scu
import spikeinterface.extractors as se 
import spikeinterface.preprocessing as spre
import spikeinterface.postprocessing as spost
import spikeinterface.sorters as ss
import spikeinterface.widgets as sw
import sys 

from probeinterface import generate_multi_columns_probe
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes
from tqdm.auto import tqdm

sys.path.append('src') 

from src.importrhdutilities import load_file

sorted_channels = {
    'D12_6': [6, 19, 26, 57, 79],
    'D13_4': [3, 4, 11, 12, 13, 15, 16, 18, 19, 23, 26, 27, 28, 29],
    'D13_8': [6, 17, 28, 75],
    'D14_4': [0, 79],
    'D14_6': [0, 32, 33, 43, 50, 51, 52, 56, 57],
}

n_s_per_min = 60
n_ms_per_s = 1000

channel_indices = np.array([
    [ 0,  1,  2,  3,  7,  6,  5,  4],
    [ 8,  9, 10, 11, 15, 14, 13, 12],
    [16, 17, 18, 19, 23, 22, 21, 20],
    [24, 25, 26, 27, 31, 30, 29, 28],
    [32, 33, 34, 35, 39, 38, 37, 36],
    [40, 41, 42, 43, 47, 46, 45, 44],
    [48, 49, 50, 51, 55, 54, 53, 52],
    [56, 57, 58, 59, 63, 62, 61, 60],
    [64, 65, 66, 67, 71, 70, 69, 68],
    [72, 73, 74, 75, 79, 78, 77, 76],
])
shank_locations = np.array([
    (115 * 0, 125 * 0), 
    (115 * 1, 125 * 1), 
    (115 * 2, 125 * 2), 
    (115 * 3, 125 * 3), 
    (115 * 4, 125 * 4), 
    (115 * 4 + 150, 125 * 4), 
    (115 * 5 + 150, 125 * 3), 
    (115 * 6 + 150, 125 * 2), 
    (115 * 7 + 150, 125 * 1),  
    (115 * 8 + 150, 125 * 0),
])

probe = generate_multi_columns_probe(
    num_columns=1, 
    num_contact_per_column=1,  
    xpitch=35,
    ypitch=35,
    contact_shapes='circle', 
    contact_shape_params={'radius': 12.5},
)
probe.set_device_channel_indices([0])


window_ms, bin_ms = 100, 1.5
isi_threshold_ms = 1.5
ms_before, ms_after = 2, 2
sorter_parameters = {
    'detect_sign': -1,
    'adjacency_radius': -1, 
    'freq_min': 300, 
    'freq_max': 3000,
    'filter': True,
    'whiten': True,  
    'clip_size': 50,
    'num_workers': 8,
    'detect_interval': 10, # default 10
}

A_active_channel_names = [
    f'A-0{channel_index:02d}' 
        for channel_index in [23, 7, 22, 6, 21, 5, 20, 4, 19, 3, 18, 2, 17, 1, 16, 0, 15, 31, 14, 30, 13, 29, 12, 28, 11, 27, 10, 26, 9, 25, 8, 24]
]
B_active_channel_names = [
    f'B-0{channel_index:02d}' 
        for channel_index in [0, 15, 1, 14, 2, 13, 3, 12] + [4, 11, 5, 10, 6, 9, 7, 8]
]
C_active_channel_names = [
    f'C-0{channel_index:02d}' 
        for channel_index in [17, 46, 19, 44, 21, 42, 23, 40, 25, 38, 27, 36, 29, 34, 31, 32, 33, 30, 35, 28, 37, 26, 39, 24, 41, 22, 43, 20, 45, 18, 47, 16]
]

active_channel_names = A_active_channel_names + B_active_channel_names + C_active_channel_names

def get_args():
    parser = argparse.ArgumentParser(description='Run Parameters')
    parser.add_argument(
        '--subject',
        type=str,
        help = 'subject name to sort',
    )
    parser.add_argument(
        '--sortdate',
        type=str,
        default=datetime.datetime.today().strftime('%Y%m%d'),
        help = 'directory to save all result',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        help = 'sorting detect threshold',
    )

    args = parser.parse_args()
    return args


def compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms):
    bins = np.arange(0, window_ms, bin_ms)
    isi = np.diff(spike_train_ms)
    if (len(isi) == 0) or (isi.min() > window_ms):
        return [], [], 0
    else:
        ys, bin_edges = np.histogram(isi, bins=bins, density=True)
        xs = bin_edges[:-1]
        rate = (isi < isi_threshold_ms).sum() / len(isi)
        return xs, ys, rate

def remove_outlier(waveforms, quantile=0.05, percentage_threshold=0.75):
    if len(waveforms) == 0:
        return waveforms
    if len(waveforms.shape) == 1:
        waveforms = waveforms[None, :]
    n_sample = waveforms.shape[1]
    bounds = np.quantile(waveforms, [quantile, 1-quantile], axis=0)
    bool_in_range = ((((waveforms >= bounds[0:1]) & (waveforms <= bounds[1:2])).sum(1) / n_sample) > percentage_threshold)
    return waveforms[bool_in_range]

def main(args):
    folder_root = f'data/processed/{args.subject}/{args.sortdate}'
    os.makedirs(folder_root, exist_ok=True)
    sorter_parameters['detect_threshold'] = args.threshold

    palette = {
        f'{args.subject}-01-before-pbs': 'orange',
        f'{args.subject}-02-after-pbs': 'green',
        f'{args.subject}-03-before-cocaine': 'blue',
        f'{args.subject}-04-after-cocaine': 'red',
    }

    recordings_folder = f'{folder_root}/recordings'
    session_info_file = f'{folder_root}/session_info.csv'
    if not os.path.isfile(session_info_file):
        session_info = []
        file_index, file_start = 0, 0
        for segment_path in (pbar := tqdm(sorted(glob.glob(f'data/raw/{args.sortdate}*/**/{args.subject}*')))):
            segment_name = segment_path.split('/')[-1]
            if segment_name in ['extra', 'injection']: continue
            pbar.set_description(segment_name)

            recording_paths = sorted(glob.glob(f'{segment_path}/*.rhd'))
            for recording_path in recording_paths:
                raw_data, data_present = load_file(recording_path)
                if data_present:
                    recording_channel_names = [channel_info['native_channel_name'] for channel_info in raw_data['amplifier_channels']]
                    sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
                    active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in active_channel_names]

                    traces = raw_data['amplifier_data'][active_channel_indices]
                    file_duration = traces.shape[1]
                    session_info.append({
                        'file_path': recording_path,
                        'file_start': file_start,
                        'file_duration': file_duration,
                        'condition': segment_name,
                        'sampling_frequency': sampling_frequency,
                    })
                    file_start += file_duration

                    recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
                    recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
                    recording = spre.common_reference(recording, reference='global', operator='median')

                    recording.save(folder=f'{recordings_folder}/file{file_index}')
                    file_index += 1

        session_info = pd.json_normalize(session_info)
        session_info.to_csv(session_info_file, index=False)
        
    session_info = pd.read_csv(session_info_file)
    n_file = len(session_info['file_path'].unique())
    recordings = [sc.load_extractor(f'{recordings_folder}/file{file_index}') for file_index in range(n_file)]
    recording = sc.concatenate_recordings(recordings)
    n_frames_per_ms = int(recording.sampling_frequency / n_ms_per_s)

            
    sortings_folder = f'{folder_root}/sortings{args.threshold}'
    waveforms_folder = f'{folder_root}/waveforms{args.threshold}'
    units_folder = f'{folder_root}/units{args.threshold}'

    for channel_index in (pbar := tqdm(range(len(recording.channel_ids)))):
        pbar.set_description(f'sorting ch{channel_index}')
        channel_sortings_folder = f'{sortings_folder}/ch{channel_index}'
        channel_waveforms_folder = f'{waveforms_folder}/ch{channel_index}'

        channel_recordings = [file_recording.channel_slice(channel_ids=file_recording.channel_ids[channel_index:channel_index+1]).set_probe(probe) for file_recording in recordings]
        channel_recording = recording.channel_slice(channel_ids=recording.channel_ids[channel_index:channel_index+1]).set_probe(probe)
        if not os.path.isfile(f'{channel_sortings_folder}/sorter_output/firings.npz'):
            ss.run_sorter(
                sorter_name='mountainsort4',
                recording=channel_recording,
                output_folder=channel_sortings_folder,
                remove_existing_folder=True,
                with_output=False,
                **sorter_parameters
            )

        channel_sorting = se.NpzSortingExtractor(f'{channel_sortings_folder}/sorter_output/firings.npz')
        # spikeinterface https://github.com/SpikeInterface/spikeinterface/pull/1378
        channel_sorting = scu.remove_excess_spikes(channel_sorting, channel_recording)
        channel_sortings = sc.split_sorting(channel_sorting, channel_recordings)
        channel_sortings = [sc.select_segment_sorting(channel_sortings, segment_indices=file_index) for file_index in range(n_file)]

        pbar.set_description(f'waveform ch{channel_index}')
        for file_index in range(n_file):
            file_waveform_folder = f'{channel_waveforms_folder}/file{file_index}'
            if not os.path.isfile(f'{channel_sortings_folder}/sorter_output/file{file_index}.npz'):
                se.NpzSortingExtractor.write_sorting(channel_sortings[file_index], f'{channel_sortings_folder}/sorter_output/file{file_index}.npz')
            if not os.path.isfile(f'{file_waveform_folder}/templates_average.npy'):
                sc.extract_waveforms(
                    channel_recordings[file_index], channel_sortings[file_index], 
                    folder=file_waveform_folder,
                    ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
                    return_scaled=False,
                    overwrite=True,
                    use_relative_path=True,
                )
        channel_waveform_extractors = [
            sc.load_waveforms(folder=f'{channel_waveforms_folder}/file{file_index}', with_recording=False, sorting=channel_sortings[file_index]) 
            for file_index in range(n_file)
        ]

        for file_index in range(n_file):
            channel_waveform_extractors[file_index].set_recording(channel_recordings[file_index])

        channel_units_folder = f'{units_folder}/ch{channel_index}'
        os.makedirs(channel_units_folder, exist_ok=True)
        subplot_size = 5
        n_plot_types = 3

        for unit_id in channel_sorting.unit_ids:
            unit_plot_file = f'{channel_units_folder}/unit{unit_id}.png'
            pbar.set_description(f'plot {args.subject} {args.threshold} unit {unit_id}')
            if True: # not os.path.isfile(unit_plot_file):
                plt.figure(figsize=(n_plot_types * subplot_size, n_file * subplot_size))
                for file_index in range(n_file):
                    condition = session_info.iloc[file_index]['condition']
                    waveforms = remove_outlier(np.squeeze(channel_waveform_extractors[file_index].get_waveforms(unit_id=unit_id)))
                    ax = plt.subplot(n_file, n_plot_types, file_index * n_plot_types + 1)
                    ax.plot(waveforms.T, color=palette[condition])
                    ax.set_title(f'[{file_index}] {waveforms.shape[0]} spikes')

                    ax = plt.subplot(n_file, n_plot_types, file_index * n_plot_types + 2)
                    if len(waveforms) > 0:
                        template = waveforms.mean(0)
                        ax.plot(template, color=palette[condition])
                    ax.set_title(f'unit {unit_id} {condition}')

                    ax = plt.subplot(n_file, n_plot_types, file_index * n_plot_types + 3)
                    spike_train_ms = channel_sortings[file_index].get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms
                    xs, ys, rate = compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms)
                    ax.bar(x=xs, height=ys, width=bin_ms, color=palette[condition], align="edge")
                    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
                    ax.set_xlabel('time (ms)')
                plt.tight_layout()
                plt.savefig(unit_plot_file)
                plt.close()
        
if __name__ == '__main__':
    args = get_args()
    main(args)