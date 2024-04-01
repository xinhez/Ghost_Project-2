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
from src.multiregion_80pin_channels import *
from src.plot import plot_unit, plot_traces

def read_recording(recording_paths, shank):
    traces = { region: [] for region in active_channel_names.keys() }
    files = []
    recording_start = 0
    for recording_path in recording_paths:
        raw_data, data_present = load_file(recording_path)
        if data_present:
            recording_channel_names = [channel['native_channel_name'] for channel in raw_data['amplifier_channels']]
            sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
            raw_traces = raw_data['amplifier_data']

            for region, region_active_channel_names in active_channel_names.items():
                region_active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in region_active_channel_names]
                region_traces = raw_traces[region_active_channel_indices]
                traces[region].append(region_traces)
            
            files.append({
                'recording_path': recording_path,
                'sampling_frequency': sampling_frequency,
                'recording_start': recording_start,  
                'recording_length': raw_traces.shape[1],
            })
            recording_start += raw_traces.shape[1] 

    recordings = {}
    for region, region_traces in traces.items():
        region_traces = np.hstack(region_traces)
        region_recording = se.NumpyRecording(traces_list=region_traces.T, sampling_frequency=sampling_frequency)
        region_recording = spre.bandpass_filter(region_recording, freq_min=300, freq_max=3000)
        region_recording = spre.common_reference(region_recording, reference='global', operator='median')
        if shank >= 0:
            region_recording = region_recording.channel_slice(channel_ids=channel_indices[shank])
        recordings[region] = region_recording

    files = pd.json_normalize(files)
    recording_duration = recording_start / sampling_frequency / n_s_per_min
    return recording_duration, recording_start, recordings, files


def sort(args, output_root, segment_paths, sorter_parameters):
    session_info_file = f'{output_root}/session_info.csv'

    recording_folder = '{output_root}/{{region}}/recording/segment{{segment}}'.format(output_root=output_root)
    if not os.path.isfile(session_info_file):
        recordings = { region: [] for region in active_channel_names.keys() }   
        session_info = []
        segment_start = 0
        segment_index = 0
        for segment_path in segment_paths:

            recording_paths = sorted(glob.glob(f'{segment_path}/*.rhd'))
            segment_duration, segment_length, segment_recordings, segment_files = read_recording(recording_paths, args.shank)
            if segment_duration >= args.min_duration:
                for region, segment_region_recording in segment_recordings.items():
                    recordings[region].append(segment_region_recording)
                    segment_region_recording.save(folder=recording_folder.format(region=region, segment=segment_index))
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

    for region in active_channel_names.keys():
        recordings = [
            sc.load_extractor(recording_folder.format(region=region, segment=segment)).set_probe(probe, in_place=True) 
            for segment in range(n_segment)
        ]
        recording = sc.concatenate_recordings(recordings).set_probe(probe, in_place=True)
        print(f'\t...Preprocessed at {recording_folder}...')

        traces_folder = f'{output_root}/{region}/traces'
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
                    title=f'{args.subject} -> {region} -> {"all" if args.shank < 0 else f"shank{args.shank}"}', 
                    savepath=trace_plot_file
                )
        print(f'\t...Plotted at {traces_folder}...')

        sorter_parameters['detect_interval'] = int(recording.sampling_frequency / n_ms_per_s * 0.3)
        sorting_folder = f'{output_root}/{region}/sorting{sorter_parameters["detect_threshold"]}'
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

        waveform_folder = f'{output_root}/{region}/waveform{sorter_parameters["detect_threshold"]}'
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
        units_folder = f'{output_root}/{region}/units{args.threshold}'
        os.makedirs(units_folder, exist_ok=True)
        for unit_id in tqdm(sorting.unit_ids):
            unit_plot_file = f'{units_folder}/{unit_id}.png'
            if not os.path.isfile(unit_plot_file):
                plot_unit(unit_id, session_info, init_date, waveform_extractors, channel_indices if args.shank < 0 else channel_indices[args.shank], savepath=unit_plot_file)
        print(f'\t...Units plotted at {units_folder}...')