import datetime
import glob 
import numpy as np
import os 
import pandas as pd
import spikeinterface.core as sc 
import spikeinterface.curation as scu
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre 
import spikeinterface.postprocessing as spost 
import spikeinterface.sorters as ss
import sys 

from tqdm.auto import tqdm

sys.path.append('src')

from src.importrhdutilities import load_file
from src.facts import *
from src.singleregion_80pin_channels import *
from src.plot import plot_unit, plot_traces

def read_recording(recording_paths, shank):
    traces = []
    files = []
    recording_start = 0
    for recording_path in recording_paths:
        raw_data, data_present = load_file(recording_path)
        if data_present:
            recording_channel_names = [channel['native_channel_name'] for channel in raw_data['amplifier_channels']]
            sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
            raw_traces = raw_data['amplifier_data']

            active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in active_channel_names]
            traces.append(raw_traces[active_channel_indices])
            
            files.append({
                'recording_path': recording_path,
                'sampling_frequency': sampling_frequency,
                'recording_start': recording_start,  
                'recording_length': raw_traces.shape[1],
            })
            recording_start += raw_traces.shape[1] 

    traces = np.hstack(traces)
    recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
    recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
    recording = spre.common_reference(recording, reference='global', operator='median')
    if shank >= 0:
        recording = recording.channel_slice(channel_ids=channel_indices[shank])

    files = pd.json_normalize(files)
    recording_duration = recording_start / sampling_frequency / n_s_per_min
    return recording_duration, recording_start, recording, files


def sort(args, output_root, segment_paths, sorter_parameters):
    session_info_file = f'{output_root}/session_info.csv'

    recording_folder = f'{output_root}/recording'
    
    if not os.path.isfile(session_info_file):
        recordings = []
        session_info = []
        segment_start = 0
        segment_index = 0
        for segment_path in segment_paths:
            recording_paths = sorted(glob.glob(f'{segment_path}/*.rhd'))
            segment_duration, segment_length, segment_recording, segment_files = read_recording(recording_paths, args.shank)
            if segment_duration >= args.min_duration:
                recordings.append(segment_recording)
                segment_recording.save(folder=f'{recording_folder}/segment{segment_index}')
                segment_files['segment_path'] = segment_path
                segment_files['segment_start'] = segment_start 
                segment_files['segment_length'] = segment_length
                session_info.append(segment_files)
                segment_start += segment_length

                segment_index += 1
        session_info = pd.concat(session_info, ignore_index=True)
        session_info.to_csv(session_info_file, index=False)
    else:
        session_info = pd.read_csv(session_info_file)

    n_segment = len(session_info['segment_path'].unique())
    if args.shank < 0:
        probe = create_multi_shank_probe(f'{output_root}/probe.png')
    else:
        probe = create_single_shank_probe(args.shank, f'{output_root}/probe.png')

    recordings = [
        sc.load_extractor(f'{recording_folder}/segment{segment}').set_probe(probe, in_place=True) 
        for segment in range(n_segment)
    ]
    recording = sc.concatenate_recordings(recordings).set_probe(probe, in_place=True)
    print(f'\t...Preprocessed at {recording_folder}...')

    traces_folder = f'{output_root}/traces'
    os.makedirs(traces_folder, exist_ok=True)
    traces = recording.get_traces().T
    n_min_plotted = 10
    duration = int(np.ceil(traces.shape[1] / recording.sampling_frequency / n_s_per_min))
    for start_min in tqdm(range(0, duration, n_min_plotted)):
        end_min = min(start_min + n_min_plotted, duration)
        trace_plot_file = f'{traces_folder}/{start_min:03d}-{end_min:03d}.png'
        if not os.path.isfile(trace_plot_file):
            plotted_traces = traces[:, int(start_min * n_s_per_min * recording.sampling_frequency) : int(end_min * n_s_per_min * recording.sampling_frequency)]
            plot_traces(
                plotted_traces, recording.sampling_frequency, 
                channel_indices if args.shank < 0 else channel_indices[args.shank:args.shank+1], 
                title=f'{args.subject} -> {"all" if args.shank < 0 else f"shank{args.shank}"}', 
                savepath=trace_plot_file
            )
    print(f'\t...Plotted at {traces_folder}...')

    sorter_parameters['detect_interval'] = int(recording.sampling_frequency / n_ms_per_s * 0.3)
    sorting_folder = f'{output_root}/sorting{sorter_parameters["detect_threshold"]}'
    if not os.path.isfile(f'{sorting_folder}/sorter_output/firings.npz'):
        ss.run_sorter(
            sorter_name='mountainsort4',
            recording=recording,
            output_folder = sorting_folder,
            remove_existing_folder=True,
            with_output=False,
            **sorter_parameters,
        )
    sorting = se.NpzSortingExtractor(f'{sorting_folder}{os.sep}sorter_output{os.sep}firings.npz')
    # spikeinterface https://github.com/SpikeInterface/spikeinterface/pull/1378
    sorting = scu.remove_excess_spikes(sorting, recording)  
    sorting = sc.split_sorting(sorting, recordings)
    sortings = [sc.select_segment_sorting(sorting, segment_indices=segment) for segment in range(n_segment)]
    for segment, segment_sorting in enumerate(sortings):
        segment_sorting_file = f'{sorting_folder}{os.sep}sorter_output{os.sep}segment{segment}_firings.npz'
        if not os.path.isfile(segment_sorting_file):
            se.NpzSortingExtractor.write_sorting(segment_sorting, segment_sorting_file)
    print(f'\t...Sorted at {sorting_folder}...')

    waveform_folder = f'{output_root}/waveform{sorter_parameters["detect_threshold"]}'
    for segment in range(n_segment):
        segment_waveform_folder = f'{waveform_folder}/segment{segment}'
        if not os.path.isfile(f'{segment_waveform_folder}{os.sep}templates_average.npy'):
            sc.extract_waveforms(
                recordings[segment], sortings[segment], 
                folder=segment_waveform_folder,
                ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
                return_scaled=False,
                overwrite=True,
                use_relative_path=True,
            )

    waveform_extractors = [sc.load_waveforms(folder=f'{waveform_folder}/segment{segment}', with_recording=False, sorting=sortings[segment]) for segment in range(n_segment)]

    for segment, recording in enumerate(recordings):
        recording.set_probe(probe, in_place=True)
        waveform_extractors[segment].set_recording(recording)
        spost.compute_unit_locations(waveform_extractors[segment], load_if_exists=False)
    print(f'\t...Waveform extracted at {waveform_folder}...')

    init_date = datetime.datetime.strptime(surgery_dates[args.subject], '%Y%m%d')
    units_folder = f'{output_root}/units{args.threshold}'
    os.makedirs(units_folder, exist_ok=True)
    for unit_id in tqdm(sorting.unit_ids):
        unit_plot_file = f'{units_folder}/{unit_id}.png'
        if not os.path.isfile(unit_plot_file):
            plot_unit(unit_id, session_info, init_date, waveform_extractors, channel_indices if args.shank < 0 else channel_indices[args.shank], savepath=unit_plot_file)
    print(f'\t...Units plotted at {units_folder}...')