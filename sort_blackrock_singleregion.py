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

sys.path.append('src')

from src.brpylib import NsxFile
from src.facts import *
from src.singleregion_32pin_channels import *
from src.plot import plot_unit, plot_traces

def read_recording(segment_path, shank):
    nsxfile = NsxFile(segment_path)
    raw_data = nsxfile.getdata()
    nsxfile.close()
    traces = np.hstack(raw_data['data'])[:n_channel] / 4
    sampling_frequency = raw_data['samp_per_s']
    
    recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
    recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
    recording = spre.common_reference(recording, reference='global', operator='median')
    if shank >= 0:
        recording = recording.channel_slice(channel_ids=channel_indices[shank])

    recording_duration = traces.shape[1] / sampling_frequency / n_s_per_min
    return recording_duration, traces.shape[1], recording


def sort(args, output_root, segment_paths, sorter_parameters):
    session_info_file = f'{output_root}/session_info.csv'

    recording_folder = f'{output_root}/recording'
    
    if not os.path.isfile(session_info_file):
        recordings = []
        session_info = []
        segment_start = 0
        segment_index = 0
        for segment_path in segment_paths:
            segment_duration, segment_length, segment_recording = read_recording(segment_path, args.shank)
            if segment_duration >= args.min_duration:
                recordings.append(segment_recording)
                if not os.path.isfile(f'{recording_folder}/segment{segment_index}/binary.json'):
                    segment_recording.save(folder=f'{recording_folder}/segment{segment_index}')
                session_info.append({
                    'segment_path': segment_path,
                    'segment_start': segment_start,
                    'segment_length': segment_length,
                })
                segment_start += segment_length
                segment_index += 1
        session_info = pd.json_normalize(session_info)
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

    traces_folder = f'{output_root}/traces'
    os.makedirs(traces_folder, exist_ok=True)
    for segment in range(n_segment):
        trace_plot_file = f'{traces_folder}/segment{segment}.png'
        segment_recording = sc.load_extractor(f'{recording_folder}/segment{segment}')
        if not os.path.isfile(trace_plot_file):
            plot_traces(
                segment_recording.get_traces().T, args.shank, segment_recording.sampling_frequency, 
                channel_indices if args.shank < 0 else channel_indices[args.shank:args.shank+1], 
                title=f'{args.subject} -> {"all" if args.shank < 0 else f"shank{args.shank}"}', 
                savepath=trace_plot_file,
            )
    print(f'\t...Plotted at {traces_folder}...')

    recordings = [
        sc.load_extractor(f'{recording_folder}/segment{segment}').set_probe(probe, in_place=True) 
        for segment in range(n_segment)
    ]
    recording = sc.concatenate_recordings(recordings).set_probe(probe, in_place=True)
    print(f'\t...Preprocessed at {recording_folder}...')

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
    for unit_id in sorting.unit_ids:
        unit_plot_file = f'{units_folder}/{unit_id}.png'
        if not os.path.isfile(unit_plot_file):
            plot_unit(unit_id, session_info, init_date, waveform_extractors, channel_indices if args.shank < 0 else channel_indices[args.shank], savepath=unit_plot_file)
    print(f'\t...Units plotted at {units_folder}...')