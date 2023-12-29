import anndata as ad
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import scanpy
import spikeinterface.preprocessing as spre
import spikeinterface.widgets as sw

from matplotlib.gridspec import GridSpec
from probeinterface import generate_multi_columns_probe
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes

# from brpylib import NsxFile

memory_limit = '20G'
n_s_per_min = 60
n_ms_per_s = 1000

blackrock_channel_indices = np.array([ 
    [ 5,  3,  1,  7,  9, 11], 
    [17, 15, 13, 19, 21, 23], 
    [29, 27, 25, 28, 26, 24], 
    [18, 20, 22, 16, 14, 12], 
    [ 6,  8, 10,  4,  2,  0],
])

intan_channel_indices = np.array([
    [24, 23, 22, 27, 26, 25], 
    [ 0, 29, 28,  3,  2,  1], 
    [ 6,  5,  4,  9,  8,  7], 
    [12, 11, 10, 15, 14, 13], 
    [18, 17, 16, 21, 20, 19],
])

def get_channels_from_the_same_shank(channel, channel_indices=blackrock_channel_indices):
    for shank in blackrock_channel_indices:
        if channel in shank:
            return shank 
    raise Exception(f'channel {channel} does not exist')

def create_probe(channel_indices):
    n_shank = len(channel_indices)
    n_channel = channel_indices.size
    shank_locations = np.array([[0, 0], [150, 200], [300, 400], [450, 200], [600, 0]])

    fig = plt.figure(figsize=(20, 10))
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
    return multi_shank_probe, fig

def preprocess_recording(recording, steps=['bp', 'cmr']):
    sampling_frequency = int(recording.sampling_frequency)
    fns = {
        'bsA': lambda recording: spre.blank_staturation(recording, abs_threshold=250, direction='both', fill_value=0),
        'bsQ': lambda recording: spre.blank_staturation(recording, quantile_threshold=0.001, direction='both', fill_value=0, chunk_size=sampling_frequency * n_s_per_min),
        'bp': lambda recording: spre.bandpass_filter(recording, freq_min=300, freq_max=3000),
        'clip': lambda recording: spre.clip(recording, a_min=-250, a_max=250),
        'cmr': lambda recording: spre.common_reference(recording, reference='global', operator='average'),
    }
    for step in steps:
        recording = fns[step](recording)
    return recording


sorter_parameters = {
    'detect_sign': -1,
    'adjacency_radius': 120, 
    'freq_min': None, 
    'freq_max': None,
    'filter': False,
    'whiten': True,  
    'clip_size': 50,
    'detect_threshold': 4,
    'detect_interval': 3, # 0.3ms 
}

ms_before = 1
ms_after = 2

def sample_indices(objects, max_n=10000):
    if len(objects) < max_n:
        return objects
    else:
        return objects[np.random.choice(np.arange(len(objects)), max_n, replace=False)]

def remove_outlier(waveforms, quantile=0.05, percentage_threshold=0.75):
    n_sample = waveforms.shape[1]
    bounds = np.quantile(waveforms, [quantile, 1-quantile], axis=0)
    bool_in_range = ((((waveforms >= bounds[0:1]) & (waveforms <= bounds[1:2])).sum(1) / n_sample) > percentage_threshold)
    return waveforms[bool_in_range]


def plot_unit(waveform_extractor, extremum_channels, sorting, unit_id, savepath, sessions=None, n_frames_per_ms=None):
    plt.rcParams.update({'font.size': 15})
    n_rows = 5
    n_cols = 3
    plt.figure(figsize=(20, 20))
    unit_extremum_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channels[unit_id]]
    unit_extremum_template = waveform_extractor.get_template(unit_id)[:, extremum_channels[unit_id]]
    deoutlier_extremum_waveforms = remove_outlier(waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channels[unit_id]])

    ax = plt.subplot(n_rows, n_cols, 1)
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:'black'})

    ax = plt.subplot(n_rows, n_cols, 2)
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

    ax = plt.subplot(n_rows, n_cols, 4)
    sw.plot_autocorrelograms(sorting, window_ms=150.0, bin_ms=1.0, unit_ids=[unit_id], axes=[ax])

    ax = plt.subplot(n_rows, n_cols, 5)
    sw.plot_isi_distribution(sorting, window_ms=150.0, bin_ms=1.0, unit_ids=[unit_id], axes=[ax])

    ax = plt.subplot(n_rows, n_cols, 7)
    ax.plot(sample_indices(unit_extremum_waveforms).T, label=unit_id, lw=0.5)
    ax.set_title(f'{unit_id} waveforms at ch {extremum_channels[unit_id]}')

    ax = plt.subplot(n_rows, n_cols, 8)
    ax.plot(unit_extremum_template.T, label=unit_id)
    ax.set_title(f'{unit_id} template at ch {extremum_channels[unit_id]}')

    ax = plt.subplot(n_rows, n_cols, 10)
    if len(deoutlier_extremum_waveforms) > 0:
        ax.plot(sample_indices(deoutlier_extremum_waveforms).T, label=unit_id, lw=0.5)
        ax.set_title(f'de-outlier {unit_id} waveforms at ch {extremum_channels[unit_id]}')
    else:
        ax.set_title(f'no waveforms left from de-outlier')

    ax = plt.subplot(n_rows, n_cols, 11)
    if len(deoutlier_extremum_waveforms) > 0:
        deoutlier_extremum_template = deoutlier_extremum_waveforms.mean(0, keepdims=True)
        ax.plot(deoutlier_extremum_template.T, label=unit_id)
        ax.set_title(f'de-outlier {unit_id} template at ch {extremum_channels[unit_id]}')
    else:
        ax.set_title(f'no template left from de-outlier')


    if sessions is not None: 
        mouse = sessions['mouse'].unique().item()
        unit_spike_train = sorting.get_unit_spike_train(unit_id)
        unit_shank_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, get_channels_from_the_same_shank(extremum_channels[unit_id])]
        unit_shank_waveforms = unit_shank_waveforms.transpose(0, 2, 1).reshape(-1, (ms_before + ms_after) * n_frames_per_ms * n_channel_per_shank)
        dates = sorted(sessions['date'])

        # Distribute spike trains to session dates.
        dated_spike_trains = {}
        for session_i in range(len(sessions)):
            session = sessions.iloc[session_i:session_i+1]
            spike_train = dated_spike_trains.get(session['date'].item(), [])
            spike_train.append(np.where(
                ( unit_spike_train >= ( session['session_start'].item() + ms_before * n_frames_per_ms ) ) &
                ( unit_spike_train <= ( session['session_start'].item() + session['session_length'].item() - ms_after * n_frames_per_ms ) ) 
            )[0])
            dated_spike_trains[session['date'].item()] = spike_train
        for date, spike_trains in dated_spike_trains.items():
            dated_spike_trains[date] = np.concatenate(spike_trains)

        # Create anndata for waveforms from the extremum shank (horizontally concatenated), labeled with the session date.
        adatas = []
        for date, spike_train in dated_spike_trains.items():
            adata = ad.AnnData(unit_shank_waveforms[sample_indices(spike_train, 1000)])
            adata.obs['date'] = date 
            adata.obs['mouse'] = mouse 
            adatas.append(adata)
        unit_shank_waveform_adata = ad.concat(adatas, index_unique='#')

        # Create anndata for waveforms from the extremum channel, labeled with the session date.
        adatas = []
        for date, spike_train in dated_spike_trains.items():
            adata = ad.AnnData(unit_extremum_waveforms[sample_indices(spike_train, 1000)])
            adata.obs['date'] = date 
            adata.obs['mouse'] = mouse 
            adatas.append(adata)
        unit_extremum_waveform_adata = ad.concat(adatas, index_unique='#')

        ax = plt.subplot(n_rows, 2, 9)
        for date_i, date in enumerate(dates):
            date_shank_waveform_adata = unit_shank_waveform_adata[unit_shank_waveform_adata.obs['date'] == date]
            ax.plot(date_shank_waveform_adata.X.mean(0))
            
        ax = plt.subplot(1, n_cols, 3)
        trace_gap = -100
        means = []
        for date_i, date in enumerate(dates):
            date_extremum_waveform_adata = unit_extremum_waveform_adata[unit_extremum_waveform_adata.obs['date'] == date]
            ax.plot(date_extremum_waveform_adata.X.mean(0) + date_i * trace_gap)

        plt.savefig(f'{savepath}0.png', bbox_inches='tight')
        plt.close()

        plt.figure(figsize=(20, 20))
        scanpy.pp.pca(unit_shank_waveform_adata, n_comps=min(50, unit_shank_waveform_adata.shape[0]-1))
        scanpy.pp.neighbors(unit_shank_waveform_adata)
        scanpy.tl.umap(unit_shank_waveform_adata)
        ax = plt.subplot(projection='3d')
        ax.set_box_aspect((1, 1, 10))
        ax.grid(False)
        shank_gap = -10
        means = []
        for date_i, date in enumerate(dates):
            date_shank_waveform_adata = unit_shank_waveform_adata[unit_shank_waveform_adata.obs['date'] == date]
            ax.scatter3D(date_shank_waveform_adata.obsm['X_umap'][:, 0], date_shank_waveform_adata.obsm['X_umap'][:, 1], date_i * shank_gap, s=2)
            means.append([date_shank_waveform_adata.obsm['X_umap'][:, 0].mean(), date_shank_waveform_adata.obsm['X_umap'][:, 1].mean(), date_i * shank_gap])
        means = np.array(means)
        ax.plot3D(means[:, 0], means[:, 1], means[:, 2])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks(np.arange(len(sessions)) * shank_gap)
        ax.set_zticklabels(dates)
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.tick_params(axis='z', pad=250, labelsize=15)
        plt.savefig(f'{savepath}1.png', bbox_inches='tight')
        plt.close()

        fig = plt.figure(figsize=(20, 20), layout="constrained")
        left = plt.imread(f'{savepath}0.png')
        right = plt.imread(f'{savepath}1.png')
        right = right[:, right.shape[1]//3+100:-right.shape[1]//4-77]

        gs = GridSpec(1, 4, figure=fig)
        ax = fig.add_subplot(gs[:, :3])
        ax.imshow(left)
        ax.set_axis_off()
        ax = fig.add_subplot(gs[:, 3:])
        ax.imshow(right)
        ax.set_axis_off()
        os.remove(f'{savepath}0.png')
        os.remove(f'{savepath}1.png')

    plt.savefig(savepath, bbox_inches='tight')
    plt.close()


########## ########## ########## ########## ########## ########## ########## ##########  
    
def plot_traces(recording, channel_indices, title, savepath, trace_gap=250, shank_gap=500, fontsize=25):
    traces = recording.get_traces().T 
    sampling_frequency = int(recording.sampling_frequency)
    n_shank, n_channel_per_shank = channel_indices.shape
    n_channel = channel_indices.size
    plt.rcParams.update({'font.size': fontsize})
    duration = traces.shape[1] / sampling_frequency / n_s_per_min
    plt.figure(figsize=(10 * duration, 50))
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
    plt.ylabel(r'200 $\mu$V gap between traces')
    plt.savefig(savepath, bbox_inches='tight')
    plt.close()

########## ########## ########## ########## ########## ########## ########## ##########  

