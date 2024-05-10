import anndata as ad 
import argparse
import glob
import matplotlib.pyplot as plt
import numpy as np 
import os
import pandas as pd 
import scanpy
import spikeinterface.core as sc 
import spikeinterface.curation as scu
import spikeinterface.extractors as se
import spikeinterface.preprocessing as spre
import spikeinterface.postprocessing as spost
import spikeinterface.sorters as ss
import spikeinterface.widgets as sw
import sys 

from probeinterface import generate_multi_columns_probe
from probeinterface.utils import combine_probes
from tqdm.auto import tqdm

sys.path.append('src')

from src.importrhdutilities import load_file

sorter_parameters = {
    'detect_sign': -1,
    'adjacency_radius': -1, 
    'freq_min': 300, 
    'freq_max': 3000,
    'filter': False,
    'whiten': True,  
    'clip_size': 50,
    'num_workers': 8,
    'detect_interval': 9, # 0.3ms
}
ms_before, ms_after = 2, 2
window_ms, bin_ms, isi_threshold_ms = 100, 1.5, 1.5
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
probes = []
for shank_channel_indices, shank_location in zip(channel_indices, shank_locations):
    probe = generate_multi_columns_probe(
        num_columns=2, num_contact_per_column=4, 
        xpitch=35, ypitch=35,
        contact_shapes='circle', contact_shape_params={'radius': 12.5}
    )
    probe.move(shank_location)
    probe.set_device_channel_indices(shank_channel_indices)

    probes.append(probe)

probe = combine_probes(probes)
probe.set_device_channel_indices(channel_indices.flatten())


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
        '--threshold',
        type=float,
        default=4.0,
        help = 'sorting detect threshold',
    )
    parser.add_argument(
        '--subject',
        type=str,
        help = 'subject to sort',
    )

    args = parser.parse_args()
    return args

def plot_autocorrelogram(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    sw.plot_autocorrelograms(waveform_extractor.sorting, window_ms=window_ms, bin_ms=bin_ms, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:color})

def plot_location(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    sw.plot_unit_locations(waveform_extractor, unit_ids=[unit_id], ax=ax)

def plot_probe_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

def plot_template_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:color})

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
    
def plot_ISI(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    n_frames_per_ms = int(waveform_extractor.sorting.sampling_frequency / n_ms_per_s)
    spike_train_ms = waveform_extractor.sorting.get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms
    xs, ys, rate = compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms)
    ax.bar(x=xs, height=ys, width=bin_ms, color=color, align="edge")
    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
    ax.set_xlabel('time (ms)')

def plot_template_extremum(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    unit_extremum_template = waveform_extractor.get_template(unit_id)[:, extremum_channel]
    ax.plot(unit_extremum_template.T, label=unit_id)
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_template(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    ax.plot(templates, label=unit_id, color=color)
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_waveforms(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    ax.plot(waveforms.T, label=unit_id, lw=0.5, color=color)
    ax.set_title(f'{unit_id} {len(waveforms) / waveform_extractor.get_total_duration():0.2f}Hz')

def plot_UMAP(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, color):
    if len(adata) > 0:
        ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], color=color, s=100)

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

def plot_unit(unit_id, waveform_extractors, file_indices, colors, savepath,
              plot_types=['autocorrelogram', 'location', 'probe_map', 'template_map', 'template_extremum', 'template', 'waveforms', 'ISI', 'UMAP'],
              subplot_size=5, min_spikes=8
              ):
    plt.rcParams.update({'font.size': 20})

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
    n_rows = len(file_indices)
    n_cols = len(plot_types)
    plt.figure(figsize=(n_cols * subplot_size, n_rows * subplot_size))
    waveforms, templates, extremum_channels, dates = [], [], [], []
    for segment in range(n_rows):
        waveform_extractor = waveform_extractors[segment]
        extremum_channel = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg')[unit_id]
        extremum_channels.append(extremum_channel)

        segment_waveforms = waveform_extractor.get_waveforms(unit_id)
        segment_templates = waveform_extractor.get_template(unit_id)
        waveforms.append(segment_waveforms.transpose(0, 2, 1).reshape(segment_waveforms.shape[0], segment_waveforms.shape[1]*segment_waveforms.shape[2]))
        templates.append(segment_templates.T.flatten())

    adata = ad.AnnData(np.vstack(waveforms))
    adata.obs['file_index'] = np.hstack([[file_indices[segment]] * len(segment_waveforms) for segment, segment_waveforms in enumerate(waveforms)])

    if len(adata) < min_spikes:
        return

    scanpy.pp.neighbors(adata, use_rep='X')
    scanpy.tl.umap(adata)

    for segment, file_index in enumerate(file_indices):
        segment_adata = adata[adata.obs['file_index'] == file_index]
        if len(segment_adata) == 0:
            continue
        for plot_i, plot_type in enumerate(plot_types):
            ax = plt.subplot(n_rows, n_cols, plot_i + segment * n_cols+1)
            plot_fns[plot_type](ax, unit_id, waveform_extractors[segment], extremum_channels[segment], waveforms[segment], templates[segment], segment_adata, colors[file_index])
            if plot_type == 'UMAP':
                ax.set_xlim(adata.obsm['X_umap'][:, 0].min(), adata.obsm['X_umap'][:, 0].max())
                ax.set_ylim(adata.obsm['X_umap'][:, 1].min(), adata.obsm['X_umap'][:, 1].max())
                ax.set_title(f'[{segment}]' + ax.get_title(), fontsize=25)
            if plot_type == 'autocorrelogram':
                ax.set_title(f'[{file_index}] unit (' + ax.get_title() + f')', fontsize=25)
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

def main(args):
    sortdate = '240506'
    subject = args.subject
    threshold = args.threshold
    sorter_parameters['detect_threshold'] = threshold

    output_root = f'data/processed/{subject}/{sortdate}'
    os.makedirs(output_root, exist_ok=True)

    print(f'{"*"*20} Processing {subject} {"*"*20}')

    recordings_folder = f'{output_root}/recordings'
    recording_paths = sorted(glob.glob(f'data/raw/{sortdate}_saline/{subject}*/{subject}-0*/*.rhd'))

    files_file = f'{output_root}/files.csv'
    if not os.path.isfile(files_file):
        files = []
        file_index = 0
        for recording_path in recording_paths:
            print(recording_path.replace('240506_saline', '240506_saline-traces-curated').replace('.rhd', '.png'))
            curated_trace_path = recording_path.replace('240506_saline', '240506_saline/240506_saline-traces-curated').replace('.rhd', '.png') 
            curated_trace_path = '/'.join(curated_trace_path.split('/')[:-2]+curated_trace_path.split('/')[-1:])
            if not os.path.isfile(curated_trace_path): continue
            condition = '-'.join(recording_path.split('/')[-2].split('-')[-2:])

            raw_data, data_present = load_file(recording_path)
            if data_present:
                recording_channel_names = [channel_info['native_channel_name'] for channel_info in raw_data['amplifier_channels']]
                sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
                active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in active_channel_names]
                traces = raw_data['amplifier_data'][active_channel_indices]

                recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
                recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
                recording = spre.common_reference(recording, reference='global', operator='median')
                recording.save(folder=f'{recordings_folder}/file{file_index}')

                files.append({
                    'file_index': file_index,
                    'file_path': recording_path,
                    'file_duration': traces.shape[1],
                    'condition': condition,
                })

                file_index += 1
        files = pd.json_normalize(files)
        files.to_csv(files_file, index=False)
    
    files = pd.read_csv(files_file)
    n_file = len(files)
    recordings = [sc.load_extractor(f'{recordings_folder}/file{file_index}').set_probe(probe) for file_index in range(n_file)]
    recording = sc.concatenate_recordings(recordings).set_probe(probe)

    sortings_folder = f'{output_root}/sortings-{threshold}'
    if not os.path.isfile(f'{sortings_folder}/sorter_output/firings.npz'):
        sorting = ss.run_sorter(
            sorter_name='mountainsort4',
            recording=recording,
            output_folder=sortings_folder,
            remove_existing_folder=True,
            with_output=True,
            **sorter_parameters
        )

    sorting = se.NpzSortingExtractor(f'{sortings_folder}/sorter_output/firings.npz')
    # spikeinterface https://github.com/SpikeInterface/spikeinterface/pull/1378
    sorting = scu.remove_excess_spikes(sorting, recording)
    sortings = sc.split_sorting(sorting, recordings)
    sortings = [sc.select_segment_sorting(sortings, segment_indices=segment) for segment in range(len(recordings))]

    waveforms_folder = f'{output_root}/waveforms-{threshold}'        
    for file_index in files['file_index'].to_list():
        file_waveform_folder = f'{waveforms_folder}/file{file_index}'
        if not os.path.isfile(f'{file_waveform_folder}/templates_average.npy'):
            sc.extract_waveforms(
                recordings[file_index], sortings[file_index], 
                folder=file_waveform_folder,
                ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
                return_scaled=False,
                overwrite=True,
                use_relative_path=True,
            )
    waveform_extractors = [sc.load_waveforms(folder=f'{waveforms_folder}/file{file_index}', with_recording=True, sorting=sortings[file_index]) for file_index in files['file_index'].to_list()]

    for waveform_extractor in waveform_extractors:
        spost.compute_unit_locations(waveform_extractor, load_if_exists=False)
    
    units_folder = f'{output_root}/units-{threshold}'
    os.makedirs(units_folder, exist_ok=True)
    for unit_id in sorting.unit_ids:
        unit_plot_file = f'{units_folder}/{unit_id}.png'
        if not os.path.isfile(unit_plot_file):
            plot_unit(unit_id, waveform_extractors, files['file_index'].to_list(), ['orange' if condition == 'before-pbs' else 'green' for condition in files['condition'].to_list()], savepath=unit_plot_file)

if __name__ == '__main__':
    args = get_args()
    main(args)