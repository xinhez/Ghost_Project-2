import matplotlib
matplotlib.use('Agg')

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

if __name__ == '__main__':
    os.makedirs('sorter_tmp', exist_ok=True)
    input_root = f'data{os.sep}raw{os.sep}MultiRegion'
    output_root = f'data{os.sep}processed{os.sep}MultiRegion5'
    sorter_parameters = {
        'detect_sign': -1,
        'adjacency_radius': -1, 
        'freq_min': 300, 
        'freq_max': 3000,
        'filter': True,
        'whiten': True,  
        'clip_size': 50,
        'detect_threshold': 5,
        'detect_interval': 10, # 0.3ms 
        'tempdir': 'sorter_tmp',
    }
    ms_before, ms_after = 1, 2
    mice = ['M9_4', 'M9_7', 'M10_1']#, 'M10_5', 'M10_8', 'M11_4']
    today = '240225'

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
        plt.show()
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

    print(active_channel_names)

    for mouse in (pbar := tqdm(mice)): 
        session_paths = glob.glob(f'{input_root}{os.sep}{mouse}{os.sep}{today}{os.sep}**')
        for session_path in session_paths:
            pbar.set_description(session_path)
            session = session_path.split(os.sep)[-1]
            session_output_folder = f'{output_root}{os.sep}{mouse}{os.sep}{session}'
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

                    if not os.path.isfile(f'{region_output_folder}{os.sep}processed{os.sep}traces.png'):
                        plot_traces(recording_processed.get_traces().T, recording_processed.sampling_frequency, channel_indices, title=f'{mouse} -> {session}', savepath=f'{region_output_folder}{os.sep}processed{os.sep}traces.png', trace_gap=150)

                recording_processed = sc.load_extractor(f'{region_output_folder}{os.sep}processed')
                
                if not os.path.isfile(f'{region_output_folder}{os.sep}sorting{os.sep}sorter_output{os.sep}firings.npz'):
                    ss.run_sorter(
                        sorter_name='mountainsort4',
                        recording=recording_processed,
                        output_folder = f'{region_output_folder}{os.sep}sorting',
                        remove_existing_folder=True,
                        with_output=True,
                        **sorter_parameters,
                    )
                sorting = se.NpzSortingExtractor(f'{region_output_folder}{os.sep}sorting{os.sep}sorter_output{os.sep}firings.npz')

                if not os.path.isfile(f'{region_output_folder}{os.sep}waveforms{os.sep}templates_average.npy'):
                    sc.extract_waveforms(
                        recording_processed, sorting, 
                        folder=f'{region_output_folder}{os.sep}waveforms',
                        ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
                        return_scaled=False,
                        overwrite=True,
                        use_relative_path=True,
                    )
                waveform_extractor = sc.load_waveforms(
                    folder=f'{region_output_folder}{os.sep}waveforms', with_recording=True, sorting=sorting
                )
                extremum_channels = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg') 
                
                units_folder = f'{region_output_folder}{os.sep}units'
                os.makedirs(units_folder, exist_ok=True)
                for unit_id in tqdm(sorting.unit_ids):
                    unit_plot_file = f'{units_folder}{os.sep}{unit_id}.png'
                    if not os.path.isfile(unit_plot_file):
                        plot_unit(waveform_extractor, extremum_channels, sorting, unit_id, channel_indices, initdate='20240128', savepath=unit_plot_file)

            session_info = pd.json_normalize(session_info)
            session_info.to_csv(f'{session_output_folder}{os.sep}session_info.csv', index=False)


