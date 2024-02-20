import matplotlib
matplotlib.use('Agg')

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
import spikeinterface.extractors as se
import spikeinterface.curation as scu
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

from src.brpylib import NsxFile


# Constants 
n_frames_per_ms = None
n_ms_per_s = 1000
n_s_per_min = 60
n_day_per_month= 30
total_month = 6
window_ms = 150
bin_ms = 2.5
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
        'detect_threshold': 5,
        'detect_interval': 3, # 0.3ms 
    },
    'mountainsort5': {

    }
}
ms_before, ms_after = 2, 2
surgery_dates = {
    '1_2': '20230720',
    '1_4': '20230704',
    '1_5': '20230627',
    '1_6': '20230630',
    '4_4': '20230804',
    '4_6': '20230804',
    '5_4': '20230807',
    '5_6': '20230807',
    '5_7': '20230805',
    '6_2': '20231019',
    '6_3': '20231016',
    '6_6': '20231016',
    '6_7': '20231016',
    '7_2': '20231019',
    '7_3': '20231027',
    '8_1': '20231126',
    '8_5': '20231126',
    '8_6': '20231126',
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
    recordings = []
    session_info = []
    session_start = 0
    for path in paths:
        nsxfile = NsxFile(path)
        raw_data = nsxfile.getdata()
        nsxfile.close()
        traces = np.hstack(raw_data['data'])[:n_channel]
        sampling_frequency = raw_data['samp_per_s'] # assume the sampling frequency is always the same for the same file format.
        duration = sum([header['data_time_s'] for header in raw_data['data_headers']])
        if duration < min_duration:
            continue 
        recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
        recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
        recording = spre.common_reference(recording, reference='global', operator='median')
        recordings.append(recording)
        session_info.append({
            'path': path,
            'sampling_frequency': sampling_frequency,
            'duration': duration,
            'session_start': session_start,
            'session_length': traces.shape[1],
        })
        session_start += traces.shape[1]
    session_info = pd.json_normalize(session_info)
    return recordings, session_info


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

def plot_autocorrelogram(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_autocorrelograms(waveform_extractor.sorting, window_ms=window_ms, bin_ms=bin_ms, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:plt.cm.turbo(lapse / total_month)})

def plot_location(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_unit_locations(waveform_extractor, unit_ids=[unit_id], ax=ax)

def plot_probe_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

def plot_template_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:plt.cm.turbo(lapse / total_month)})

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
    spike_train_ms = waveform_extractor.sorting.get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms
    xs, ys, rate = compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms, isi_threshold_ms)
    ax.bar(x=xs, height=ys, width=bin_ms, color=plt.cm.turbo(lapse / total_month), align="edge")
    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
    ax.set_xlabel('time (ms)')

def plot_template(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    ax.plot(templates, label=unit_id, color=plt.cm.turbo(lapse / total_month))
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_waveforms(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    ax.plot(waveforms.T, label=unit_id, lw=0.5, color=plt.cm.turbo(lapse / total_month))
    ax.set_title(f'{unit_id} waveforms at ch {extremum_channel}')

def plot_UMAP(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, lapse):
    if len(adata) > 0:
        ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], color=plt.cm.turbo(lapse / total_month), s=20)
        ax.set_title(f'{lapse} months')

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
    
def plot_unit(
        unit_id, session_info, init_date, waveform_extractors, savepath,
        plot_types=['autocorrelogram', 'location', 'probe_map', 'template_map', 'template', 'waveforms', 'ISI', 'UMAP'],
        subplot_size=5
    ):

    plot_fns = {
        'autocorrelogram': plot_autocorrelogram, 
        'location': plot_location, 
        'probe_map': plot_probe_map, 
        'template_map': plot_template_map, 
        'ISI': plot_ISI, 
        'template': plot_template, 
        'waveforms': plot_waveforms, 
        'UMAP': plot_UMAP,
    }

    plt.figure(figsize=(len(plot_types) * subplot_size, len(session_info) * subplot_size))
    waveforms, templates, extremum_channels, dates = [], [], [], []
    for segment, session_path in enumerate(session_info['path']):
        waveform_extractor = waveform_extractors[segment]
        extremum_channel = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg')[unit_id]
        extremum_shank = get_shank(extremum_channel)
        extremum_channels.append(extremum_channel)

        segment_waveforms = sample_objects(waveform_extractor.get_waveforms(unit_id)[:, :, channel_indices[extremum_shank]], max_n=100)
        segment_templates = waveform_extractor.get_template(unit_id)[:, channel_indices[extremum_shank]]
        waveforms.append(segment_waveforms.transpose(0, 2, 1).reshape(segment_waveforms.shape[0], segment_waveforms.shape[1]*segment_waveforms.shape[2]))
        templates.append(segment_templates.T.flatten())

        match = re.search(r'(\d{4})(\d{2})(\d{2})', session_path)
        dates.append(datetime.datetime.strptime(match.group(0), '%Y%m%d'))

    adata = ad.AnnData(np.vstack(waveforms))
    adata.obs['segment'] = np.hstack([[segment] * len(segment_waveforms) for segment, segment_waveforms in enumerate(waveforms)])
    scanpy.pp.neighbors(adata, use_rep='X')
    scanpy.tl.umap(adata)

    for segment, segment_date in enumerate(dates):
        segment_adata = adata[adata.obs['segment'] == segment]
        if len(segment_adata) == 0:
            continue
        lapse = round(((segment_date - init_date).days) / n_day_per_month)
        for plot_i, plot_type in enumerate(plot_types):
            ax = plt.subplot(len(session_info), len(plot_types), plot_i + segment * len(plot_types)+1)
            plot_fns[plot_type](ax, unit_id, waveform_extractors[segment], extremum_channels[segment], waveforms[segment], templates[segment], segment_adata, lapse)
            if plot_type == 'UMAP':
                ax.set_xlim(adata.obsm['X_umap'][:, 0].min(), adata.obsm['X_umap'][:, 0].max())
                ax.set_ylim(adata.obsm['X_umap'][:, 1].min(), adata.obsm['X_umap'][:, 1].max())
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

def main(args):
    folder_root = f'{args.folder}/{args.subject}/{args.savedate}'
    print(f'Saving results to {folder_root}...')

    recording_folder = f'{folder_root}/recording'
    if not os.path.isfile(f'{recording_folder}/binary.json'):
        recording_paths = sorted(glob.glob(f'{args.data}/{args.subject}/**/*.{args.nsx}'))
        recordings, session_info = read_recordings(recording_paths, args.min_duration, channel_indices.size)   
        recording = sc.concatenate_recordings(recordings)
        recording.save(folder=recording_folder)
        session_info.to_csv(f'{folder_root}/session_info.csv', index=False)
    else:
        session_info = pd.read_csv(f'{folder_root}/session_info.csv')
        recording_paths = session_info['path'].tolist()
        recordings, _ = read_recordings(recording_paths, args.min_duration, channel_indices.size)   

    probe = create_probe(channel_indices, shank_locations, f'{recording_folder}/probe.png')
    recording = sc.load_extractor(recording_folder).set_probe(probe, in_place=True)
    global n_frames_per_ms
    n_frames_per_ms = recording.sampling_frequency / n_ms_per_s
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
    
    sorting = sc.split_sorting(sorting, recordings)
    sortings = [sc.select_segment_sorting(sorting, segment_indices=segment) for segment in range(len(session_info))]
    print('\t...Sorted...')

    waveform_folder = f'{folder_root}/waveform'
    for segment in range(len(session_info)):
        segment_waveform_folder = f'{waveform_folder}/{segment}'
        if not os.path.isfile(f'{segment_waveform_folder}{os.sep}templates_average.npy'):
            recordings[segment].set_probe(probe, in_place=True).save(folder='tmp-recording')
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

    waveform_extractors = [sc.load_waveforms(folder=f'{waveform_folder}/{segment}', with_recording=False, sorting=sortings[segment]) for segment in range(len(session_info))]

    for segment, recording in enumerate(recordings):
        recording.set_probe(probe, in_place=True)
        waveform_extractors[segment].set_recording(recording)
        spost.compute_unit_locations(waveform_extractors[segment], load_if_exists=False)
    print('\t...Waveform extracted...')

    init_date = datetime.datetime.strptime(surgery_dates[args.subject], '%Y%m%d')
    units_folder = f'{folder_root}/units'
    os.makedirs(units_folder, exist_ok=True)
    for unit_id in tqdm(sorting.unit_ids):
        unit_plot_file = f'{units_folder}/{unit_id}.png'
        if not os.path.isfile(unit_plot_file):
            plot_unit(unit_id, session_info, init_date, waveform_extractors, savepath=unit_plot_file)

if __name__ == '__main__':
    args = get_args()
    main(args)