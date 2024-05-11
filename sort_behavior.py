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
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes
from tqdm.auto import tqdm 

sys.path.append('.')

from src.importrhdutilities import load_file

n_s_per_min = 60 
n_ms_per_s = 1000
n_total_trial = 100

channels_by_region = {
    'CA1': ['A-023', 'A-007', 'A-022', 'A-006', 'A-021', 'A-005', 'A-020', 'A-004', 'A-019', 'A-003', 'A-018', 'A-002', 'A-017', 'A-001', 'A-016', 'A-000', 'A-015', 'A-031', 'A-014', 'A-030', 'A-013', 'A-029', 'A-012', 'A-028', 'A-011', 'A-027', 'A-010', 'A-026', 'A-009', 'A-025', 'A-008', 'A-024', 'B-000', 'B-015', 'B-001', 'B-014', 'B-002', 'B-013', 'B-003', 'B-012'],
    'M1': ['B-004', 'B-011', 'B-005', 'B-010', 'B-006', 'B-009', 'B-007', 'B-008', 'C-017', 'C-046', 'C-019', 'C-044', 'C-021', 'C-042', 'C-023', 'C-040', 'C-025', 'C-038', 'C-027', 'C-036', 'C-029', 'C-034', 'C-031', 'C-032', 'C-033', 'C-030', 'C-035', 'C-028', 'C-037', 'C-026', 'C-039', 'C-024', 'C-041', 'C-022', 'C-043', 'C-020', 'C-045', 'C-018', 'C-047', 'C-016'],
}   

channel_indices = np.array([
    [ 0,  1,  2,  3,  7, 6, 5, 4], 
    [ 8,  9, 10, 11, 15, 14, 13, 12], 
    [16, 17, 18, 19, 23, 22, 21, 20], 
    [24, 25, 26, 27, 31, 30, 29, 28],
    [32, 33, 34, 35, 39, 38, 37, 36], 
])
shank_locations = [(0, 300), (100, 150), (200, 0), (300, 150), (400, 300)]

window_ms, bin_ms = 100, 1.5
ms_before, ms_after = 2, 2
sorter_parameters = {
    'detect_sign': -1,
    'adjacency_radius': -1, 
    'freq_min': 300, 
    'freq_max': 3000,
    'filter': False,
    'whiten': True,  
    'clip_size': 50,
    'detect_interval': 9, # 0.3ms 
    'num_workers': 28,
}

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


def plot_autocorrelogram(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    sw.plot_autocorrelograms(waveform_extractor.sorting, window_ms=window_ms, bin_ms=bin_ms, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:plt.cm.tab10(segment % 10)})

def plot_location(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    sw.plot_unit_locations(waveform_extractor, unit_ids=[unit_id], ax=ax)

def plot_probe_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

def plot_template_map(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:plt.cm.tab10(segment % 10)})

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
    
def plot_ISI(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment, isi_threshold_ms=1.5):
    n_frames_per_ms = int(waveform_extractor.sorting.sampling_frequency / n_ms_per_s)
    spike_train_ms = waveform_extractor.sorting.get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms
    xs, ys, rate = compute_isi_violation_rate(spike_train_ms, window_ms, bin_ms, isi_threshold_ms)
    ax.bar(x=xs, height=ys, width=bin_ms, color=plt.cm.tab10(segment % 10), align="edge")
    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
    ax.set_xlabel('time (ms)')

def plot_template_extremum(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    unit_extremum_template = waveform_extractor.get_template(unit_id)[:, extremum_channel]
    ax.plot(unit_extremum_template.T, label=unit_id)
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_template(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    ax.plot(templates, label=unit_id, color=plt.cm.tab10(segment % 10))
    ax.set_title(f'{unit_id} template at ch {extremum_channel}')

def plot_waveforms(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    ax.plot(waveforms.T, label=unit_id, lw=0.5, color=plt.cm.tab10(segment % 10))
    ax.set_title(f'{unit_id} waveforms at ch {extremum_channel}')

def plot_UMAP(ax, unit_id, waveform_extractor, extremum_channel, waveforms, templates, adata, segment):
    if len(adata) > 0:
        ax.scatter(adata.obsm['X_umap'][:, 0], adata.obsm['X_umap'][:, 1], color=plt.cm.tab10(segment % 10), s=100)
        ax.set_title(f'[{segment}]', fontsize=25)

def get_shank(channel_id, channel_indices):
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

def plot_unit(unit_id, waveform_extractors, channel_indices, savepath,
              plot_types=['autocorrelogram', 'location', 'probe_map', 'template_map', 'template_extremum', 'template', 'waveforms', 'ISI', 'UMAP'],
              subplot_size=5, min_spikes=10):
    plt.rcParams.update({'font.size': 25})
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
    n_rows = len(waveform_extractors)
    n_cols = len(plot_types)
    plt.figure(figsize=(n_cols * subplot_size, n_rows * subplot_size))
    waveforms, templates, extremum_channels = [], [], []
    for segment_index in range(len(waveform_extractors)):
        waveform_extractor = waveform_extractors[segment_index]
        extremum_channel = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg')[unit_id]
        extremum_channel = np.where(waveform_extractors[segment_index].channel_ids == extremum_channel)[0].item()
        extremum_shank = get_shank(extremum_channel, channel_indices)
        extremum_channel_indices = channel_indices[extremum_shank]
        extremum_channels.append(extremum_channel)

        segment_waveforms = sample_objects(waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channel_indices], max_n=100)
        segment_templates = waveform_extractor.get_template(unit_id)[:, extremum_channel_indices]
        waveforms.append(segment_waveforms.transpose(0, 2, 1).reshape(segment_waveforms.shape[0], segment_waveforms.shape[1]*segment_waveforms.shape[2]))
        templates.append(segment_templates.T.flatten())

    adata = ad.AnnData(np.vstack(waveforms))
    adata.obs['segment'] = np.hstack([[segment] * len(segment_waveforms) for segment, segment_waveforms in enumerate(waveforms)])

    if len(adata) < min_spikes:
        return

    scanpy.pp.neighbors(adata, use_rep='X')
    scanpy.tl.umap(adata)

    for segment_index in range(len(waveform_extractors)):
        segment_adata = adata[adata.obs['segment'] == segment_index]
        if len(segment_adata) == 0:
            continue
        for plot_i, plot_type in enumerate(plot_types):
            ax = plt.subplot(n_rows, n_cols, plot_i + segment_index * n_cols+1)
            plot_fns[plot_type](ax, unit_id, waveform_extractors[segment_index], extremum_channels[segment_index], waveforms[segment_index], templates[segment_index], segment_adata, segment_index)
            if plot_type == 'UMAP':
                ax.set_xlim(adata.obsm['X_umap'][:, 0].min(), adata.obsm['X_umap'][:, 0].max())
                ax.set_ylim(adata.obsm['X_umap'][:, 1].min(), adata.obsm['X_umap'][:, 1].max())
    plt.tight_layout()
    plt.savefig(savepath)
    plt.close()

def get_args():
    parser = argparse.ArgumentParser(description='Run Parameters')
    parser.add_argument(
        '--subject',
        type=str,
        help='subject name to sort',
    )
    parser.add_argument(
        '--threshold',
        type=float,
        default=5.5,
        help='sorting detect threshold',
    )
    args = parser.parse_args()
    return args

def main(args):
    sorter_parameters['detect_threshold'] = args.threshold

    output_root = f'data/processed/{args.subject}'
    print('*'*20, output_root, args.threshold, '*'*20)

    session_info_file = f'{output_root}/info.csv'

    if not os.path.isfile(session_info_file):
        session_info = []
        segment_start = 0
        segment_paths = sorted(glob.glob(f'data/raw/intan/{args.subject}/{args.subject}*'))
        for segment_index, segment_path in (pbar := tqdm(enumerate(segment_paths))):
            pbar.set_description(f'Reading segment {segment_index + 1}/{len(segment_paths)}')
            recording_paths = sorted(glob.glob(f'{segment_path}/*.rhd'))

            segment_traces = { 'CA1': [], 'M1': [] }
            triggers = []
            for recording_path in recording_paths:
                raw_data, data_present = load_file(recording_path)
                if data_present:
                    recorded_channels = [channel_info['native_channel_name'] for channel_info in raw_data['amplifier_channels']]
                    for region, region_channels in channels_by_region.items():
                        segment_traces[region].append(
                            raw_data['amplifier_data'][
                                [recorded_channels.index(region_channel) for region_channel in region_channels]
                        ])
                    triggers.append(raw_data['board_dig_in_data'][0])
                    sampling_frequency = int(raw_data['frequency_parameters']['amplifier_sample_rate'])
            triggers = np.hstack(triggers)
            segment_traces = { region: np.hstack(region_traces) for region, region_traces in segment_traces.items() }

            trigger_ups = np.where(np.diff(triggers.astype(int), prepend=triggers[0].astype(int)) == 1)[0]
            trigger_downs = np.where(np.diff(triggers.astype(int), prepend=triggers[0].astype(int)) == -1)[0]

            # If recording started with trigger up, skip this trigger and go to the next. 
            if trigger_ups[0] > trigger_downs[0]:
                trigger_downs = trigger_downs[1:]

            trigger_durations = np.round((trigger_downs - trigger_ups[:len(trigger_downs)]) / sampling_frequency)
            trial_triggers = np.where(trigger_durations >= 2)[0]
            assert len(trial_triggers) in [n_total_trial, n_total_trial+1]
            trial_triggers = trial_triggers[:n_total_trial]
            trial_starts = trigger_ups[trial_triggers]

            # A delay of idle_duration is set at the beginning of the experiment.
            idle_duration = 10 # s 
            # Need to look at preset_duration seconds before the first trial to get the ratser plot.
            preset_duration = 5 # s
            # Each tigger is up for signal_duraiton seconds to signal the start of a trial.
            signal_duration = 2 # s

            assert len(trial_triggers) == n_total_trial 
            assert trigger_durations[trial_triggers[0]] == idle_duration + signal_duration
            assert (len(trigger_ups) - len(trigger_durations)) <= 1

            trial_starts[0] = trial_starts[0] + idle_duration * sampling_frequency
            trial_ends = np.append(trial_starts[1:], trigger_ups[-1])

            experiment_start = trial_starts[0] - preset_duration * sampling_frequency
            experiment_end = trial_ends[-1]
            segment_duration = experiment_end - experiment_start

            trial_starts = trial_starts - experiment_start
            trial_ends = trial_ends - experiment_start

            segment_traces = { region: region_traces[:, experiment_start:experiment_end] for region, region_traces in segment_traces.items() }

            for region, region_traces in segment_traces.items():
                if not os.path.isfile(f'{output_root}/{region}/recordings/segment{segment_index}/binary.json'):
                    recording = se.NumpyRecording(traces_list=region_traces.T, sampling_frequency=sampling_frequency)
                    recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
                    recording = spre.common_reference(recording, operator='median', reference='global')
                    recording.save(folder=f'{output_root}/{region}/recordings/segment{segment_index}')

            session_info.append({
                'subject': args.subject,
                'segment_index': segment_index,
                'segment_path': segment_path,
                'segment_start': segment_start,
                'segment_duration': segment_duration,
                'sampling_frequency': sampling_frequency,
                'trial_starts': trial_starts.tolist(),
                'trial_ends': trial_ends.tolist(),
            })
            segment_start += segment_duration
        session_info = pd.json_normalize(session_info)
        session_info.to_csv(session_info_file, index=False)
    else:
        session_info = pd.read_csv(session_info_file)
        n_segment = len(session_info)
        probe = create_probe(
            channel_indices, 
            shank_locations, 
            n_rows=4, n_cols=2, 
            inter_electrode_distance=30, 
            electrode_radius=10, savepath=f'{output_root}/probe'
        )
        for region in channels_by_region.keys():
            recordings = [sc.load_extractor(f'{output_root}/{region}/recordings/segment{segment_index}').set_probe(probe) for segment_index in range(n_segment)]
            recording = sc.concatenate_recordings(recordings).set_probe(probe)
            print(recording)
            
            sortings_folder = f'{output_root}/{region}/sortings-{args.threshold}'
            if not os.path.isfile(f'{sortings_folder}/sorter_output/firings.npz'):
                ss.run_sorter(
                    sorter_name='mountainsort4',
                    recording=recording,
                    output_folder = sortings_folder,
                    remove_existing_folder=True,
                    with_output=False,
                    **sorter_parameters,
                )
            sorting = se.NpzSortingExtractor(f'{sortings_folder}/sorter_output/firings.npz')
            sorting = scu.remove_excess_spikes(sorting, recording)
            sortings = sc.split_sorting(sorting, recordings)
            sortings = [sc.select_segment_sorting(sortings, segment_indices=segment) for segment in range(len(recordings))]

            waveforms_folder = f'{output_root}/{region}/waveforms-{args.threshold}'
            for segment_index in range(n_segment):
                segment_waveform_folder = f'{waveforms_folder}/segment{segment_index}'
                if not os.path.isfile(f'{segment_waveform_folder}/templates_average.npy'):
                    sc.extract_waveforms(
                        recordings[segment_index], sortings[segment_index], 
                        folder=segment_waveform_folder,
                        ms_before=ms_before, ms_after=ms_after, max_spikes_per_unit=None,
                        return_scaled=False,
                        overwrite=True,
                        use_relative_path=True,
                    )
            waveform_extractors = [sc.load_waveforms(folder=f'{waveforms_folder}/segment{segment_index}', with_recording=False, sorting=sortings[segment_index]) for segment_index in range(n_segment)]
            for segment_index in range(n_segment):
                waveform_extractors[segment_index].set_recording(recordings[segment_index])
                spost.compute_unit_locations(waveform_extractors[segment_index], load_if_exists=False)

            units_folder = f'{output_root}/{region}/units-{args.threshold}'
            os.makedirs(units_folder, exist_ok=True)
            for unit_id in (pbar := tqdm (sorting.unit_ids)):
                pbar.set_description(f'plotting {unit_id}/{len(sorting.unit_ids)}')
                unit_plot_file = f'{units_folder}/{unit_id}.png'
                if not os.path.isfile(unit_plot_file):
                    plot_unit(unit_id, waveform_extractors, channel_indices, savepath=unit_plot_file)

if __name__ == '__main__':
    args = get_args()
    main(args)