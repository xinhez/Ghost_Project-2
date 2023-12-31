import anndata as ad
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

memory_limit = '20G'
n_s_per_min = 60
n_ms_per_s = 1000
min_recording_duration = 10

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

ms_before, ms_after = 1, 2
window_ms, bin_ms = 100, 1.5
isi_threshold_ms = 1.5

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

def sample_objects(objects, max_n=None):
    if max_n is None:
        return objects
    if len(objects) < max_n:
        return objects
    else:
        return objects[np.random.choice(np.arange(len(objects)), max_n, replace=False)]

def remove_outlier(waveforms, quantile=0.05, percentage_threshold=0.75):
    n_sample = waveforms.shape[1]
    bounds = np.quantile(waveforms, [quantile, 1-quantile], axis=0)
    bool_in_range = ((((waveforms >= bounds[0:1]) & (waveforms <= bounds[1:2])).sum(1) / n_sample) > percentage_threshold)
    return waveforms[bool_in_range]

def compute_isi_violation_rate(spike_train_ms):
    bins = np.arange(0, window_ms, bin_ms)
    isi = np.diff(spike_train_ms)
    if (len(isi) == 0) or (isi.min() > window_ms):
        return [], [], 0
    else:
        ys, bin_edges = np.histogram(isi, bins=bins, density=True)
        xs = bin_edges[:-1]
        rate = (isi < isi_threshold_ms).sum() / len(isi)
        return xs, ys, rate

def get_unit_session_spike_train(sorting, unit_id, session, as_indices):
    n_frames_per_ms = sorting.sampling_frequency / n_ms_per_s
    unit_spike_train = sorting.get_unit_spike_train(unit_id=unit_id)
    session_spike_train_indices = np.where(
        ( unit_spike_train >= ( session['session_start'].item() + ms_before * n_frames_per_ms ) ) &
        ( unit_spike_train <= ( session['session_start'].item() + session['session_length'].item() - ms_after * n_frames_per_ms ) ) 
    )[0]
    if as_indices:
        return session_spike_train_indices
    else:
        return unit_spike_train[session_spike_train_indices] - session['session_start'].item()
    
def split_unit_spike_train_indicies_by_session(sorting, unit_id, sessions):
    session_spike_trains = {}
    for session_i in range(len(sessions)):
        session = sessions[session_i: session_i+1]
        session_spike_trains[session_i] = get_unit_session_spike_train(sorting, unit_id, session, as_indices=True)
    return session_spike_trains

def get_unit_shank_waveforms(waveform_extractor, extremum_channels, channel_indices, unit_id):
    n_channel_per_shank = channel_indices.shape[1]
    n_frames_per_ms = int(waveform_extractor.sampling_frequency / n_ms_per_s)
    unit_shank_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, get_channels_from_the_same_shank(extremum_channels[unit_id], channel_indices)]
    unit_shank_waveforms = unit_shank_waveforms.transpose(0, 2, 1).reshape(-1, (ms_before + ms_after) * n_frames_per_ms * n_channel_per_shank)
    return unit_shank_waveforms

def get_session_waveforms_adata(waveforms, session_spike_trains, max_count_per_session):
    adatas = []
    for session_i, spike_train in session_spike_trains.items():
        adata = ad.AnnData(waveforms[sample_objects(spike_train, max_count_per_session)])
        adata.obs['session_i'] = session_i 
        adatas.append(adata)
    return ad.concat(adatas, index_unique='#')

def erase_3D_pane(ax):
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.pane.set_edgecolor('w')
    ax.yaxis.pane.set_edgecolor('w')
    ax.zaxis.pane.set_edgecolor('w')

def collate_lateral_figures(left, right, split=(7, 2)):
    fig = plt.figure(figsize=(20, 20), layout="constrained")
    gs = GridSpec(1, sum(split), figure=fig)
    ax = fig.add_subplot(gs[:, :split[0]])
    ax.imshow(left)
    ax.set_axis_off()
    ax = fig.add_subplot(gs[:, split[0]:])
    ax.imshow(right)
    ax.set_axis_off()

def plot_unit(waveform_extractor, extremum_channels, sorting, unit_id, channel_indices, savepath, sessions=None):
    plt.rcParams.update({'font.size': 15})
    n_rows = 5
    n_cols = 3
    plt.figure(figsize=(20, 20))
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    unit_extremum_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channels[unit_id]]
    unit_extremum_template = waveform_extractor.get_template(unit_id)[:, extremum_channels[unit_id]]
    deoutlier_extremum_waveforms = remove_outlier(waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channels[unit_id]])
    n_frames_per_ms = sorting.sampling_frequency / n_ms_per_s
    spike_train_ms = sorting.get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms

    ax = plt.subplot(n_rows, n_cols, 1)
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:'black'})

    ax = plt.subplot(n_rows, n_cols, 2)
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

    ax = plt.subplot(n_rows, n_cols, 4)
    sw.plot_autocorrelograms(sorting, window_ms=100, bin_ms=1.5, unit_ids=[unit_id], axes=[ax])
    ax.set_xlabel('time (ms)')

    ax = plt.subplot(n_rows, n_cols, 5)
    xs, ys, rate = compute_isi_violation_rate(spike_train_ms)
    ax.bar(x=xs, height=ys, width=bin_ms, color="gray", align="edge")
    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
    ax.set_xlabel('time (ms)')

    ax = plt.subplot(n_rows, n_cols, 7)
    ax.plot(sample_objects(unit_extremum_waveforms).T, label=unit_id, lw=0.5)
    ax.set_title(f'{unit_id} waveforms at ch {extremum_channels[unit_id]}')

    ax = plt.subplot(n_rows, n_cols, 8)
    ax.plot(unit_extremum_template.T, label=unit_id)
    ax.set_title(f'{unit_id} template at ch {extremum_channels[unit_id]}')

    ax = plt.subplot(n_rows, n_cols, 10)
    if len(deoutlier_extremum_waveforms) > 0:
        ax.plot(sample_objects(deoutlier_extremum_waveforms).T, label=unit_id, lw=0.5)
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
        session_spike_trains = split_unit_spike_train_indicies_by_session(sorting, unit_id, sessions)
        unit_extremum_waveform_adata = get_session_waveforms_adata(unit_extremum_waveforms, session_spike_trains, max_count_per_session=1000)
        unit_shank_waveforms = get_unit_shank_waveforms(waveform_extractor, extremum_channels, channel_indices, unit_id)
        unit_shank_waveform_adata = get_session_waveforms_adata(unit_shank_waveforms, session_spike_trains, max_count_per_session=1000)

        scanpy.pp.pca(unit_shank_waveform_adata, n_comps=min(50, unit_shank_waveform_adata.shape[0]-1))
        scanpy.pp.neighbors(unit_shank_waveform_adata)
        scanpy.tl.umap(unit_shank_waveform_adata)

        ax = plt.subplot(n_rows, 2, 9)
        for session_i in range(len(sessions)):
            session_shank_waveform_adata = unit_shank_waveform_adata[unit_shank_waveform_adata.obs['session_i'] == session_i]
            ax.plot(session_shank_waveform_adata.X.mean(0), color=plt.cm.turbo(session_i / len(sessions)))
            
        ax = plt.subplot(1, n_cols, 3)
        trace_gap = -100
        means = []
        for session_i in range(len(sessions)):
            session_extremum_waveform_adata = unit_extremum_waveform_adata[unit_extremum_waveform_adata.obs['session_i'] == session_i]
            ax.plot(session_extremum_waveform_adata.X.mean(0) + session_i * trace_gap, color=plt.cm.turbo(session_i / len(sessions)))
        ax.set_yticks(np.arange(0, trace_gap * len(sessions), trace_gap), sessions['date'])
        ax.yaxis.tick_right()

        plt.savefig(f'{savepath}0.png', bbox_inches='tight') # Hack to remove 3D plot excessive margins.
        plt.close()

        plt.figure(figsize=(20, 20))
        ax = plt.subplot(projection='3d')
        ax.set_box_aspect((1, 1, 10))
        erase_3D_pane(ax)

        shank_gap = -10
        means = []
        for session_i in range(len(sessions)):
            session_shank_waveform_adata = unit_shank_waveform_adata[unit_shank_waveform_adata.obs['session_i'] == session_i]
            ax.scatter3D(session_shank_waveform_adata.obsm['X_umap'][:, 0], session_shank_waveform_adata.obsm['X_umap'][:, 1], session_i * shank_gap, s=2, color=plt.cm.turbo(session_i / len(sessions)))
            means.append([session_shank_waveform_adata.obsm['X_umap'][:, 0].mean(), session_shank_waveform_adata.obsm['X_umap'][:, 1].mean(), session_i * shank_gap])
        means = np.array(means)
        ax.plot3D(means[:, 0], means[:, 1], means[:, 2])
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks(np.arange(len(sessions)) * shank_gap, sessions['date'])
        ax.set_xlabel('UMAP1')
        ax.set_ylabel('UMAP2')
        ax.tick_params(axis='z', pad=250, labelsize=15)
        plt.savefig(f'{savepath}1.png', bbox_inches='tight')  # Hack to remove 3D plot excessive margins.
        plt.close()

        left = plt.imread(f'{savepath}0.png')
        right = plt.imread(f'{savepath}1.png')
        right = right[:, right.shape[1]//3+100:-right.shape[1]//4-77]
        collate_lateral_figures(left, right)
        os.remove(f'{savepath}0.png')
        os.remove(f'{savepath}1.png')

    plt.savefig(savepath, bbox_inches='tight')
    plt.close()

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

def plot_isi_by_session(sorting, unit_id, sessions, savepath=None):
    plt.figure(figsize=(15, 15))
    ax = plt.subplot(projection='3d')
    ax.set_box_aspect((4, 20, 2))
    erase_3D_pane(ax)

    ax.set_xlabel('time (ms)')
    ax.set_ylabel(f'sessions ISI Violation Rate ({isi_threshold_ms}ms)', labelpad=180)
    ax.set_ylim(0, len(sessions))
    ax.set_zlabel('frequency', labelpad=20)
    ax.tick_params(axis='z', pad=10)
    ax.set_zlim(0, 0.05)

    n_frames_per_ms = sorting.sampling_frequency / n_ms_per_s

    ytick_labels = []
    for session_i in range(len(sessions)):
        session = sessions.iloc[session_i:session_i+1]
        date = session['date'].item()
        spike_train_ms = get_unit_session_spike_train(sorting, unit_id, session, as_indices=False) / n_frames_per_ms

        xs, ys, rate = compute_isi_violation_rate(spike_train_ms)
        ax.bar(xs, ys, zs=session_i, zdir='y', width=bin_ms, align="edge", label=date, color=plt.cm.turbo(session_i/len(sessions)))
        ytick_labels.append(f'{date} {rate*100:0.1f}%')

    ax.set_yticks(np.arange(len(sessions))+1, ytick_labels)
    ax.tick_params(axis='y', direction='out', pad=50, labelrotation=-15)

    ax.set_title(f'ISI unit {unit_id}')
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    plt.show()
    plt.close()