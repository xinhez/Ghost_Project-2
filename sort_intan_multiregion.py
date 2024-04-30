import matplotlib
matplotlib.use('Agg')

import argparse
import datetime
import glob
import matplotlib.pyplot  as plt
import numpy as np
import os
import pandas as pd
import spikeinterface.core as sc
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
import spikeinterface.sorters as ss
import sys 

from probeinterface import generate_multi_columns_probe, ProbeGroup
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes
from tqdm.auto import tqdm

sys.path.append('src')

from src.importrhdutilities import load_file
from utils import *


ms_before, ms_after = 2, 2

channel_indices = np.array([
    [ 0,  1,  2,  3,  7, 6, 5, 4], #4,  5,  6,  7],
    [ 8,  9, 10, 11, 15, 14, 13, 12], #12, 13, 14, 15],
    [16, 17, 18, 19, 23, 22, 21, 20], #20, 21, 22, 23],
    [24, 25, 26, 27, 31, 30, 29, 28], #28, 29, 30, 31],
    [32, 33, 34, 35, 39, 38, 37, 36], #36, 37, 38, 39],
])
shank_locations = [(0, 300), (100, 150), (200, 0), (300, 150), (400, 300)]

def create_probe(channel_indices, shank_locations, n_rows, n_cols, inter_electrode_distance, electrode_radius, savepath=None):
    plt.rcParams.update({'font.size': 8})
    n_shank = len(channel_indices)
    n_channel = channel_indices.size
    
    plt.figure(figsize=(20, 10))
    ax = plt.gca()

    probes = []
    for shank_channel_indices, shank_location in zip(channel_indices, shank_locations):
        probe = generate_multi_columns_probe(
            num_columns=n_cols, num_contact_per_column=n_rows, 
            xpitch=inter_electrode_distance, ypitch=inter_electrode_distance,
            contact_shapes='circle', contact_shape_params={'radius': electrode_radius}
        )
        probe.move(shank_location)
        probe.set_device_channel_indices(shank_channel_indices)

        plot_probe(probe, with_device_index=True, ax=ax)
        probes.append(probe)
        
    multi_shank_probe = combine_probes(probes)
    multi_shank_probe.set_device_channel_indices(channel_indices.flatten())

    plt.xlim(-100, 600)
    plt.ylim(-150, 600)
    plt.title(f'Probe - {n_channel}ch - {n_shank}shanks')
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    # plt.show()
    plt.close()
    return multi_shank_probe

A_active_channel_names = [
    f'A-0{channel_index:02d}' 
        for channel_index in [23, 7, 22, 6, 21, 5, 20, 4, 19, 3, 18, 2, 17, 1, 16, 0, 15, 31, 14, 30, 13, 29, 12, 28, 11, 27, 10, 26, 9, 25, 8, 24]
]
B1_active_channel_names = [
    f'B-0{channel_index:02d}' 
        for channel_index in [0, 15, 1, 14, 2, 13, 3, 12]
]

B2_active_channel_names = [
    f'B-0{channel_idnex:02d}'
        for channel_idnex in [4, 11, 5, 10, 6, 9, 7, 8]
]
C_active_channel_names = [
    f'C-0{channel_index:02d}' 
        for channel_index in [17, 46, 19, 44, 21, 42, 23, 40, 25, 38, 27, 36, 29, 34, 31, 32, 33, 30, 35, 28, 37, 26, 39, 24, 41, 22, 43, 20, 45, 18, 47, 16]
]


active_channel_names = {
    'CA1': A_active_channel_names + B1_active_channel_names,
    'M1': B2_active_channel_names + C_active_channel_names,
}

def get_args():
    """ Get command line arguments. """
    parser = argparse.ArgumentParser(description='Run Parameters')
    parser.add_argument(
        '--subject',
        type=str,
        help = 'subject name to sort',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        help = 'sorting detect threshold',
    )
    parser.add_argument(
        '--sortdate',
        type=str,
        default=datetime.datetime.today().strftime('%Y%m%d'),
        help = 'directory to save all result',
    )
    args = parser.parse_args()
    return args


def plot_traces(traces, sampling_frequency, channel_indices, title, savepath, session_w=10, trace_gap=250, shank_gap=500, fontsize=25):
    n_shank, n_channel_per_shank = channel_indices.shape
    n_channel = channel_indices.size
    plt.rcParams.update({'font.size': fontsize})
    duration = traces.shape[1] / sampling_frequency / n_s_per_min
    plt.figure(figsize=(session_w * duration, 50))
    plt.title(f'{title} : {duration:0.2f} min')

    for shank_i, shank in enumerate(channel_indices):
        for channel_i, channel in enumerate(shank):
            y_baseline = trace_gap * (channel_i + shank_i * n_channel_per_shank) + shank_gap * shank_i
            
            plt.plot(traces[channel]+y_baseline)
            plt.text(len(traces[channel]), y_baseline - fontsize, f'ch{channel}')

    xticks_labels = list(range(round(duration) + 1))
    xticks_locs = [min * sampling_frequency * n_s_per_min for min in xticks_labels]
    plt.ylim(-trace_gap, n_channel * trace_gap + (n_shank - 1) * shank_gap)
    plt.xticks(ticks=xticks_locs, labels=xticks_labels)
    plt.xlabel('min')
    plt.ylabel(rf'{trace_gap} $\mu$V gap between traces')
    plt.savefig(savepath, bbox_inches='tight')
    plt.close()
    
def main(args):
    os.makedirs('sorter_tmp', exist_ok=True)
    input_root = f'data{os.sep}raw{os.sep}MultiRegion'
    output_root = f'data{os.sep}processed{os.sep}MultiRegion'

    sorter_parameters = {
        'detect_sign': -1,
        'adjacency_radius': -1, 
        'freq_min': 300, 
        'freq_max': 3000,
        'filter': False,
        'whiten': True,  
        'clip_size': 50,
        'detect_threshold': args.threshold,
        'detect_interval': 9, # 0.3ms
        'num_workers': 8,
    }
        
    session_paths = glob.glob(f'{input_root}{os.sep}{args.subject}{os.sep}{args.sortdate}{os.sep}**')

    for session_path in (pbar := tqdm(session_paths)):
        pbar.set_description(session_path)
        session = session_path.split(os.sep)[-1]
        session_output_folder = f'{output_root}{os.sep}{args.subject}{os.sep}{session}'
        os.makedirs(session_output_folder, exist_ok=True)
        
        recording_paths = sorted(glob.glob(f'{session_path}{os.sep}*.rhd'))

        traces = { 'CA1': [], 'M1': [] }
        session_info = []
        session_length = 0
        for recording_path in recording_paths:
            raw_data, data_present = load_file(recording_path)
            if data_present:
                recording_channel_names = [channel_info['native_channel_name'] for channel_info in raw_data['amplifier_channels']]
                for region, region_active_channel_names in active_channel_names.items():
                    recording_active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in region_active_channel_names]
                    raw_trace = raw_data['amplifier_data'][recording_active_channel_indices]
                    traces[region].append(raw_trace)
                sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
                session_info.append({
                    'file': recording_path,
                    'file_start': session_length,
                    'file_length': raw_data['amplifier_data'].shape[1],
                    'sampling_frequency': sampling_frequency,
                    'session': session,
                })
                session_length += raw_data['amplifier_data'].shape[1] 
        session_info = pd.json_normalize(session_info)
        session_info.to_csv(f'{session_output_folder}{os.sep}session_info.csv', index=False)
            
        for region, region_traces in traces.items():
            region_output_folder = f'{session_output_folder}{os.sep}{region}'
            os.makedirs(region_output_folder, exist_ok=True)

            if not os.path.isfile(f'{region_output_folder}{os.sep}processed{os.sep}binary.json'):
                traces[region] = np.hstack(region_traces) 
                recording = se.NumpyRecording(traces_list=traces[region].T, sampling_frequency=sampling_frequency)

                multi_shank_probe = create_probe(
                    channel_indices, 
                    shank_locations, 
                    n_rows=4, n_cols=2, 
                    inter_electrode_distance=30, 
                    electrode_radius=10, savepath=f'{region_output_folder}{os.sep}probe'
                )
                recording.set_probe(multi_shank_probe, in_place=True)

                recording_processed = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
                recording_processed = spre.common_reference(recording_processed, reference='global', operator='median')
                recording_processed.save(folder=f'{region_output_folder}{os.sep}processed')

            recording_processed = sc.load_extractor(f'{region_output_folder}{os.sep}processed')

            if not os.path.isfile(f'{region_output_folder}{os.sep}traces.png'):
                plot_traces(recording_processed.get_traces().T, recording_processed.sampling_frequency, channel_indices, title=f'{args.subject} -> {session}', savepath=f'{region_output_folder}{os.sep}traces.png', trace_gap=150)
                
            traces_folder = f'{region_output_folder}{os.sep}session_traces'
            os.makedirs(traces_folder, exist_ok=True)
            for session_i in tqdm(range(len(session_info))):
                session_file = session_info.loc[session_i, 'file'].split(os.sep)[-1].replace('rhd', 'png')
                session_trace_file = f'{traces_folder}{os.sep}{session_file}'
                if not os.path.isfile(session_trace_file):
                    session_start = session_info.loc[session_i, 'file_start']
                    session_end = session_start + session_info.loc[session_i, 'file_length']
                    plot_traces(recording_processed.get_traces(start_frame=session_start, end_frame=session_end).T, recording_processed.sampling_frequency, channel_indices, title=f'{args.subject} -> {session}', savepath=session_trace_file, trace_gap=150, session_w=50)
            
            # sorting_folder = f'{region_output_folder}{os.sep}sorting{args.threshold}'
            # if not os.path.isfile(f'{sorting_folder}{os.sep}sorter_output{os.sep}firings.npz'):
            #     ss.run_sorter(
            #         sorter_name='mountainsort4',
            #         recording=recording_processed,
            #         output_folder = sorting_folder,
            #         remove_existing_folder=True,
            #         with_output=True,
            #         **sorter_parameters,
            #     )
            # sorting = se.NpzSortingExtractor(f'{sorting_folder}{os.sep}sorter_output{os.sep}firings.npz')

            # waveforms_folder = f'{region_output_folder}{os.sep}waveforms{args.threshold}'
            # if not os.path.isfile(f'{waveforms_folder}{os.sep}templates_average.npy'):
            #     sc.extract_waveforms(
            #         recording_processed, sorting, 
            #         folder=waveforms_folder,
            #         ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
            #         return_scaled=False,
            #         overwrite=True,
            #         use_relative_path=True,
            #     )
            # waveform_extractor = sc.load_waveforms(
            #     folder=f'{waveforms_folder}', with_recording=True, sorting=sorting
            # )
            # extremum_channels = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg') 
            
            # units_folder = f'{region_output_folder}{os.sep}units{args.threshold}'
            # os.makedirs(units_folder, exist_ok=True)
            # for unit_id in tqdm(sorting.unit_ids):
            #     unit_plot_file = f'{units_folder}{os.sep}{unit_id}.png'
            #     if not os.path.isfile(unit_plot_file):
            #         plot_unit(waveform_extractor, extremum_channels, sorting, unit_id, channel_indices, initdate='20240101', savepath=unit_plot_file)


if __name__ == '__main__':
    args = get_args()
    main(args)