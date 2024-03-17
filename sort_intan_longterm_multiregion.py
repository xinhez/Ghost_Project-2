import anndata as ad
import argparse
import datetime
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
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

n_day_per_week = 7 
n_s_per_min = 60
n_ms_per_s = 1000
window_ms = 100 
bin_ms = 1.5
total_week = 20
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
channel_indices = np.array([
    [ 0,  1,  2,  3,  7, 6, 5, 4], #4,  5,  6,  7],
    [ 8,  9, 10, 11, 15, 14, 13, 12], #12, 13, 14, 15],
    [16, 17, 18, 19, 23, 22, 21, 20], #20, 21, 22, 23],
    [24, 25, 26, 27, 31, 30, 29, 28], #28, 29, 30, 31],
    [32, 33, 34, 35, 39, 38, 37, 36], #36, 37, 38, 39],
])
shank_locations = np.array([(0, 300), (100, 150), (200, 0), (300, 150), (400, 300)])
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
surgery_dates = {
    'M9_4':  '20240211',
    'M9_7':  '20240211',
    'M10_6': '20240222',
    'M10_8': '20240128',
    'M11_4': '20240128',
}

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
        '--sorter',
        type=str,
        default='mountainsort4',
        help = 'sorting algorithm',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        help = 'sorting detect threshold',
    )
    parser.add_argument(
        '--min_duration',
        type=int,
        default=10,
        help = 'minimum duration (s) to consider as valid data',
    )
    args = parser.parse_args()
    return args

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

def read_recording(recording_paths):
    traces = { 'CA1':[], 'M1':[] }
    files = []
    file_start = 0
    for recording_path in recording_paths:
        raw_data, data_present = load_file(recording_path)
        if data_present:
            recording_channel_names = [channel_info['native_channel_name'] for channel_info in raw_data['amplifier_channels']]
            sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
            for region in traces:
                region_active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in active_channel_names[region]]
                traces[region].append(raw_data['amplifier_data'][region_active_channel_indices])
            file_duration = raw_data['amplifier_data'].shape[1]
            files.append({
                'file_path': recording_path,
                'file_start': file_start,
                'file_duration': file_duration,
                'sampling_frequency': sampling_frequency,
            })
            file_start += file_duration

    recording = { 'CA1':[], 'M1':[] }
    for region in traces:
        region_traces = np.hstack(traces[region])
        recording[region] = se.NumpyRecording(traces_list=region_traces.T, sampling_frequency=sampling_frequency)
        recording[region] = spre.bandpass_filter(recording[region], freq_min=300, freq_max=3000)
        recording[region] = spre.common_reference(recording[region], reference='global', operator='median')

    files = pd.json_normalize(files)
    duration = files['file_duration'].sum() / sampling_frequency / n_s_per_min
    return recording, duration, files

def plot_autocorrelogram(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_autocorrelograms(waveform_extractor.sorting, window_ms=window_ms, bin_ms=bin_ms, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:plt.cm.turbo(lapse / total_week)})

def plot_location(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_unit_locations(waveform_extractor, unit_ids=[unit_id], ax=ax)

def plot_probe_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

def plot_template_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:plt.cm.turbo(lapse / total_week)})

def compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms, isi_threshold_ms):
    bins = np.arange(0, window_ms, bin_ms)
    isi = np.diff(spike_train_ms)
    if (len(isi) == 0) or (isi.min() > window_ms):
        return [], [], 0
    else:
        ys, bin_edges = np.histogram(isi, bins=bins, density=True)
        xs = bin_edges[:-1]
        rate = (isi < isi_threshold_ms).sum() / len(isi)
        return xs, ys, rate
    
def plot_ISI(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse, isi_threshold_ms=1.5):
    n_frames_per_ms = int(waveform_extractor.sorting.sampling_frequency / n_ms_per_s)
    spike_train_ms = waveform_extractor.sorting.get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms
    xs, ys, rate = compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms, isi_threshold_ms)
    ax.bar(x=xs, height=ys, width=bin_ms, color=plt.cm.turbo(lapse / total_week), align="edge")
    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
    ax.set_xlabel('time (ms)')

def plot_template_extremum(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    unit_extremum_template = waveform_extractor.get_template(unit_id)[:, extremum_channel]
    ax.plot(unit_extremum_template.T, label=unit_id)
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_template(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    ax.plot(templates, label=unit_id, color=plt.cm.turbo(lapse / total_week))
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_waveforms(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    ax.plot(waveforms.T, label=unit_id, lw=0.5, color=plt.cm.turbo(lapse / total_week))
    ax.set_title(f'{unit_id} waveforms at ch {extremum_channel}')

def plot_UMAP(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    if len(adata) > 0:
        ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], color=plt.cm.turbo(lapse / total_week), s=100)
        ax.set_title(f'{lapse:0.2f} weeks')

def get_shank(channel_id):
    for shank, shank_channels in enumerate(channel_indices):
        if channel_id in shank_channels:
            return shank 
    raise Exception

def sample_objects(objects, max_n=None):
    if max_n is None:
        return objects
    if len(objects) < max_n:
        return objects
    else:
        return objects[np.random.choice(np.arange(len(objects)), max_n, replace=False)]

def plot_unit(unit_id, session_info, init_date, waveform_extractors, savepath,
              plot_types=['autocorrelogram', 'location', 'probe_map', 'template_map', 'template_extremum', 'template', 'waveforms', 'ISI', 'UMAP'],
              subplot_size=5, min_spikes=10
              ):

    plot_fns = {
        'autocorrelogram': plot_autocorrelogram, 
        'location': plot_location, 
        'probe_map': plot_probe_map, 
        'template_map': plot_template_map, 
        'ISI': plot_ISI, 
        'template_extremum': plot_template_extremum,
        'template': plot_template, 
        'waveforms': plot_waveforms, 
        'UMAP': plot_UMAP,
    }
    n_rows = len(session_info['segment_path'].unique())
    n_cols = len(plot_types)
    plt.figure(figsize=(n_cols * subplot_size, n_rows * subplot_size))
    waveforms, templates, extremum_channels, dates = [], [], [], []
    for segment, segment_path in enumerate(session_info['segment_path'].unique()):
        waveform_extractor = waveform_extractors[segment]
        extremum_channel = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg')[unit_id]
        extremum_shank = get_shank(extremum_channel)
        extremum_channels.append(extremum_channel)

        segment_waveforms = sample_objects(waveform_extractor.get_waveforms(unit_id)[:, :, channel_indices[extremum_shank]], max_n=100)
        segment_templates = waveform_extractor.get_template(unit_id)[:, channel_indices[extremum_shank]]
        waveforms.append(segment_waveforms.transpose(0, 2, 1).reshape(segment_waveforms.shape[0], segment_waveforms.shape[1]*segment_waveforms.shape[2]))
        templates.append(segment_templates.T.flatten())

        dates.append(datetime.datetime.strptime(segment_path.split('/')[-1].split('_')[-2], '%y%m%d'))

    adata = ad.AnnData(np.vstack(waveforms))
    adata.obs['segment'] = np.hstack([[segment] * len(segment_waveforms) for segment, segment_waveforms in enumerate(waveforms)])

    if len(adata) < min_spikes:
        return

    scanpy.pp.neighbors(adata, use_rep='X')
    scanpy.tl.umap(adata)

    for segment, segment_date in enumerate(dates):
        segment_adata = adata[adata.obs['segment'] == segment]
        if len(segment_adata) == 0:
            continue
        lapse = ((segment_date - init_date).days) / n_day_per_week
        for plot_i, plot_type in enumerate(plot_types):
            ax = plt.subplot(n_rows, n_cols, plot_i + segment * n_cols+1)
            plot_fns[plot_type](ax, unit_id, waveform_extractors[segment], extremum_channels[segment], waveforms[segment], templates[segment], segment_adata, lapse)
            if plot_type == 'UMAP':
                ax.set_xlim(adata.obsm['X_umap'][:, 0].min(), adata.obsm['X_umap'][:, 0].max())
                ax.set_ylim(adata.obsm['X_umap'][:, 1].min(), adata.obsm['X_umap'][:, 1].max())
                ax.set_title(f'[{segment}]' + ax.get_title(), fontsize=50)
            if plot_type == 'autocorrelogram':
                ax.set_title(f'unit (' + ax.get_title() + f') {round(lapse)} wk', fontsize=50)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

def main(args):
    sorter_parameters['detect_threshold'] = args.threshold

    folder_root = f'data/processed/longterm/{args.subject}/{args.sortdate}'
    os.makedirs(folder_root, exist_ok=True)
    print(f'Saving results to {folder_root}...')

    segment_paths = sorted(glob.glob(f'data/raw/MultiRegion/{args.subject}/**/**'))

    recordings = { 'CA1':[], 'M1':[] }
    if os.path.isfile(f'{folder_root}/session_info.csv'):
        session_info = pd.read_csv(f'{folder_root}/session_info.csv')
        for segment_path in sorted(session_info['segment_path'].unique()):
            recording_paths = session_info[session_info['segment_path'] == segment_path]['file_path']
            recording, duration, segment_files = read_recording(recording_paths)
            assert duration >=args.min_duration
            for region in recordings:
                recordings[region].append(recording[region])
    else:
        session_info = []
        segment_start = 0
        for segment_path in tqdm(segment_paths):
            recording_paths = sorted(glob.glob(f'{segment_path}/*.rhd'))
            recording, duration, segment_files = read_recording(recording_paths)
            if duration >= args.min_duration:
                for region in recordings:
                    recordings[region].append(recording[region])

                segment_duration = segment_files['file_duration'].sum()
                segment_files['segment_path'] = segment_path
                segment_files['segment_duration'] = segment_duration
                segment_files['segment_start'] = segment_start
                session_info.append(segment_files)
                segment_start += segment_duration
        session_info = pd.concat(session_info, ignore_index=True)
        session_info.to_csv(f'{folder_root}/session_info.csv', index=False)

    for region, region_recordings in recordings.items():
        region_folder = f'{folder_root}/{region}'

        recording_folder = f'{region_folder}/recording'
        if not os.path.isfile(f'{recording_folder}/binary.json'):
            recording = sc.concatenate_recordings(region_recordings)
            recording.save(folder=recording_folder)

        recording = sc.load_extractor(recording_folder)
        probe = create_probe(
            channel_indices, shank_locations, 
            n_rows=4, n_cols=2, inter_electrode_distance=30, electrode_radius=10, 
            savepath=f'{region_folder}{os.sep}probe'
        )
        recording.set_probe(probe, in_place=True)
        print('*'*20, f'Preprocessing at {recording_folder}')

        sorting_folder = f'{region_folder}/sorting{args.threshold}'
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
        sorting = sc.split_sorting(sorting, region_recordings)
        sortings = [sc.select_segment_sorting(sorting, segment_indices=segment) for segment in range(len(region_recordings))]
        print('*'*20, f'Sorting at {sorting_folder}')

        waveform_folder = f'{region_folder}/waveform{args.threshold}'
        for segment in range(len(region_recordings)):
            segment_waveform_folder = f'{waveform_folder}/{segment}'
            if not os.path.isfile(f'{segment_waveform_folder}/templates_average.npy'):
                region_recordings[segment].set_probe(probe, in_place=True).save(folder='tmp-recording')
                se.NpzSortingExtractor.write_sorting(sortings[segment], 'tmp-sorting')

                segment_recording = sc.load_extractor('tmp-recording')
                segment_sorting = se.NpzSortingExtractor('tmp-sorting.npz')
                sc.extract_waveforms(
                    segment_recording, segment_sorting, 
                    folder=segment_waveform_folder,
                    ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
                    return_scaled=False,
                    overwrite=True,
                    use_relative_path=True,
                )
                shutil.rmtree('tmp-recording')
                os.remove('tmp-sorting.npz')
        waveform_extractors = [
            sc.load_waveforms(folder=f'{waveform_folder}/{segment}', with_recording=False, sorting=sortings[segment]) 
            for segment in range(len(region_recordings))
        ]

        for segment, recording in enumerate(region_recordings):
            recording.set_probe(probe, in_place=True)
            waveform_extractors[segment].set_recording(recording)
            spost.compute_unit_locations(waveform_extractors[segment], load_if_exists=False)
        print('*'*20, f'Waveforms at {waveform_folder}')

        init_date = datetime.datetime.strptime(surgery_dates[args.subject], '%Y%m%d')
        units_folder = f'{region_folder}/units{args.threshold}'
        os.makedirs(units_folder, exist_ok=True)
        for unit_id in tqdm(sorting.unit_ids):
            unit_plot_file = f'{units_folder}/{unit_id}.png'
            if not os.path.isfile(unit_plot_file):
                plot_unit(unit_id, session_info, init_date, waveform_extractors, savepath=unit_plot_file)
            

if __name__ == "__main__":
    args = get_args()
    main(args)