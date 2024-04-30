import anndata as ad
import datetime
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
import spikeinterface.widgets as sw

from matplotlib.gridspec import GridSpec
from probeinterface import generate_multi_columns_probe, ProbeGroup
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes

from brpylib import NsxFile
from records import *

def read_nsx(recording_path, trace_scale=4):
    nsxfile = NsxFile(recording_path)
    raw_data = nsxfile.getdata()
    sampling_frequency = raw_data['samp_per_s']
    raw_traces = np.hstack(raw_data['data']) / trace_scale
    nsxfile.close()
    return raw_traces, sampling_frequency

def get_channels_from_the_same_shank(channel, channel_indices=blackrock_channel_indices):
    for shank in channel_indices:
        if channel in shank:
            return shank 
    raise Exception(f'channel {channel} does not exist')

def create_probe(channel_indices, savepath=None):
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
    plt.show()
    plt.close()
    return multi_shank_probe

def create_probegroup(channel_indices, savepath=None):
    plt.rcParams.update({'font.size': 10})
    n_shank = len(channel_indices)
    n_channel = channel_indices.size

    plt.figure(figsize=(20, 10))
    ax = plt.gca()

    probegroup = ProbeGroup()
    for shank_channel_indices, shank_location in zip(channel_indices, shank_locations):
        probe = generate_multi_columns_probe(
            num_columns=2, num_contact_per_column=3, 
            xpitch=50, ypitch=50,
            contact_shapes='circle', contact_shape_params={'radius': 10}
        )
        probe.move(shank_location)
        probe.set_device_channel_indices(shank_channel_indices)

        plot_probe(probe, with_device_index=True, ax=ax)
        probegroup.add_probe(probe)

    plt.xlim(-100, 700)
    plt.ylim(-150, 550)
    plt.title(f'Probe - {n_channel}ch - {n_shank}shanks')
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    plt.show()
    plt.close()
    return probegroup

def preprocess_recording(recording, steps=['bp', 'cmr']):
    sampling_frequency = int(recording.sampling_frequency)
    fns = {
        'bsA': lambda recording: spre.blank_staturation(recording, abs_threshold=250, direction='both', fill_value=0),
        'bsQ': lambda recording: spre.blank_staturation(recording, quantile_threshold=0.001, direction='both', fill_value=0, chunk_size=sampling_frequency * n_s_per_min),
        'bp': lambda recording: spre.bandpass_filter(recording, freq_min=300, freq_max=3000),
        'clip': lambda recording: spre.clip(recording, a_min=-250, a_max=250),
        'cmr': lambda recording: spre.common_reference(recording, reference='global', operator='median'),
    }
    for step in steps:
        recording = fns[step](recording)
    return recording

def read_sorted_results(sorted_folder, read_sessions, by_group=False):
    recording_processed = sc.load_extractor(f'{sorted_folder}{os.sep}processed')

    if by_group:
        sorting = se.NpzSortingExtractor(f'{sorted_folder}{os.sep}sorting-by-group{os.sep}sorter_output{os.sep}firings.npz')
    else:
        sorting = se.NpzSortingExtractor(f'{sorted_folder}{os.sep}sorting{os.sep}sorter_output{os.sep}firings.npz')
    sorting = scu.remove_excess_spikes(sorting, recording_processed) # spikeinterface https://github.com/SpikeInterface/spikeinterface/pull/1378

    if by_group:
        waveform_extractor = sc.load_waveforms(
            folder=f'{sorted_folder}{os.sep}waveforms-by-group', with_recording=True, sorting=sorting
        )
    else:
        waveform_extractor = sc.load_waveforms(
            folder=f'{sorted_folder}{os.sep}waveforms', with_recording=True, sorting=sorting
        )
    extremum_channels = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg') 
    
    if read_sessions:
        sessions = pd.read_csv(f'{sorted_folder}{os.sep}sessions.csv').sort_values(by='date')  
        return recording_processed, sorting, waveform_extractor, extremum_channels, sessions
    else:
        return recording_processed, sorting, waveform_extractor, extremum_channels

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
    for session_i in sessions.index:
        session = sessions[session_i: session_i+1]
        session_spike_trains[session_i] = get_unit_session_spike_train(sorting, unit_id, session, as_indices=True)
    return session_spike_trains

def get_unit_shank_waveforms(waveform_extractor, extremum_channels, channel_indices, unit_id):
    n_channel_per_shank = channel_indices.shape[1]
    n_frames_per_ms = int(waveform_extractor.sampling_frequency / n_ms_per_s)
    unit_shank_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, get_channels_from_the_same_shank(extremum_channels[unit_id], channel_indices)]
    unit_shank_waveforms = unit_shank_waveforms.transpose(0, 2, 1).reshape(-1, (ms_before + ms_after) * n_frames_per_ms * n_channel_per_shank)
    return unit_shank_waveforms

def get_waveforms_by_session_adata(waveforms, session_spike_trains, max_count_per_session=None):
    adatas = []
    for session_i, spike_train in session_spike_trains.items():
        if max_count_per_session is not None:
            spike_train = sample_objects(spike_train, max_count_per_session)
        adata = ad.AnnData(waveforms[spike_train])
        adata.obs['session_i'] = session_i 
        adatas.append(adata)
    return ad.concat(adatas, index_unique='#')

def get_unit_extremum_waveforms_by_session_adata(unit_id, sorting, waveform_extractor, extremum_channels, sessions, max_count_per_session=None):
    unit_extremum_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channels[unit_id]]
    spike_trains_by_session = split_unit_spike_train_indicies_by_session(sorting, unit_id, sessions)
    return get_waveforms_by_session_adata(unit_extremum_waveforms, spike_trains_by_session, max_count_per_session=max_count_per_session)

def get_unit_shank_waveforms_by_session_adata(unit_id, channel_indices, sorting, waveform_extractor, extremum_channels, sessions, max_count_per_session=None):
    spike_trains_by_session = split_unit_spike_train_indicies_by_session(sorting, unit_id, sessions)
    unit_shank_waveforms = get_unit_shank_waveforms(waveform_extractor, extremum_channels, channel_indices, unit_id)
    unit_shank_waveform_by_session_adata = get_waveforms_by_session_adata(unit_shank_waveforms, spike_trains_by_session, max_count_per_session=max_count_per_session)
    scanpy.pp.pca(unit_shank_waveform_by_session_adata, n_comps=min(50, unit_shank_waveform_by_session_adata.shape[0]-1))
    scanpy.pp.neighbors(unit_shank_waveform_by_session_adata)
    scanpy.tl.umap(unit_shank_waveform_by_session_adata)
    return unit_shank_waveform_by_session_adata

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

def plot_UMAP_by_session(sessions, unit_shank_waveform_adata, plotted_session_indices, initdate, savepath=None, colormap=plt.cm.turbo):        
    plt.figure(figsize=(20, 20))
    ax = plt.subplot(projection='3d')
    ax.set_box_aspect((1, 1, 10))
    erase_3D_pane(ax)

    shank_gap = -10
    means = []
    ztick_labels = []
    for plot_index, session_i in enumerate(plotted_session_indices):
        session_shank_waveform_adata = unit_shank_waveform_adata[unit_shank_waveform_adata.obs['session_i'] == session_i]
        
        ax.scatter3D(session_shank_waveform_adata.obsm['X_umap'][:, 0], session_shank_waveform_adata.obsm['X_umap'][:, 1], plot_index * shank_gap, s=2, color=plt.cm.turbo(session_i / len(sessions)))
        
        means.append([session_shank_waveform_adata.obsm['X_umap'][:, 0].mean(), session_shank_waveform_adata.obsm['X_umap'][:, 1].mean(), plot_index * shank_gap])
        
        session_date = sessions['date'][session_i]
        session_lapse = round(((datetime.datetime.strptime(str(session_date), '%Y%m%d') - datetime.datetime.strptime(initdate, '%Y%m%d')).days) / n_day_per_month)
        ztick_labels.append(f'[{session_i}] {session_date} ({session_lapse} mo)')

    means = np.array(means)
    ax.plot3D(means[:, 0], means[:, 1], means[:, 2])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks(np.arange(len(plotted_session_indices)) * shank_gap, ztick_labels)
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    ax.tick_params(axis='z', pad=320, labelsize=15)
    for plot_index, ticklabel in enumerate(ax.get_zticklabels()):
        ticklabel.set_color(colormap(plotted_session_indices[plot_index]/len(sessions)))
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')  # Hack to remove 3D plot excessive margins.
    plt.show()
    plt.close()

def plot_isi_by_session(unit_id, sorting, sessions, plotted_session_indices, initdate, savepath=None, colormap=plt.cm.turbo):
    plt.figure(figsize=(15, 15))
    ax = plt.subplot(projection='3d')
    ax.set_box_aspect((4, 20, 2))
    erase_3D_pane(ax)

    ax.set_xlabel('time (ms)', labelpad=10)
    ax.set_ylabel(f'sessions ISI Violation Rate ({isi_threshold_ms}ms)', labelpad=280)
    ax.set_ylim(0, len(plotted_session_indices))
    ax.set_zlabel('frequency', labelpad=30)
    ax.tick_params(axis='z', pad=10)
    ax.set_zlim(0, 0.05)

    n_frames_per_ms = sorting.sampling_frequency / n_ms_per_s

    ytick_labels = []
    for plot_index, session_i in enumerate(plotted_session_indices):
        session = sessions.iloc[session_i:session_i+1]
        session_date = session['date'].item()
        spike_train_ms = get_unit_session_spike_train(sorting, unit_id, session, as_indices=False) / n_frames_per_ms

        xs, ys, rate = compute_isi_violation_rate(spike_train_ms)
        ax.bar(xs, ys, zs=plot_index, zdir='y', width=bin_ms, align="edge", label=session_date, color=colormap(session_i/len(sessions)))

        session_lapse = round(((datetime.datetime.strptime(str(session_date), '%Y%m%d') - datetime.datetime.strptime(initdate, '%Y%m%d')).days) / n_day_per_month)

        ytick_labels.append(f'[{session_i}] {session_date} {rate*100:0.1f}% ({session_lapse} mo)')

    ax.set_yticks(np.arange(len(plotted_session_indices))+1, ytick_labels)
    ax.tick_params(axis='y', direction='out', pad=120)
    for plot_index, ticklabel in enumerate(ax.get_yticklabels()):
        ticklabel.set_color(colormap(plotted_session_indices[plot_index]/len(sessions)))

    ax.set_title(f'ISI unit {unit_id}')
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')

    plt.show()
    plt.close()

def plot_unit(waveform_extractor, extremum_channels, sorting, unit_id, channel_indices, initdate, savepath, sessions=None):
    plt.rcParams.update({'font.size': 15})
    n_rows = 5
    n_cols = 3
    plt.figure(figsize=(20, 20))
    plt.subplots_adjust(wspace=0.4, hspace=0.4)
    
    n_frames_per_ms = sorting.sampling_frequency / n_ms_per_s
    
    ax = plt.subplot(n_rows, n_cols, 1)
    sw.plot_unit_templates(waveform_extractor, unit_ids=[unit_id], axes=[ax], unit_colors={unit_id:'black'})

    ax = plt.subplot(n_rows, n_cols, 2)
    sw.plot_unit_probe_map(waveform_extractor, unit_ids=[unit_id], axes=[ax])

    ax = plt.subplot(n_rows, n_cols, 4)
    sw.plot_autocorrelograms(sorting, window_ms=100, bin_ms=1.5, unit_ids=[unit_id], axes=[ax])
    ax.set_xlabel('time (ms)')

    ax = plt.subplot(n_rows, n_cols, 5)
    spike_train_ms = sorting.get_unit_spike_train(unit_id=unit_id) / n_frames_per_ms
    xs, ys, rate = compute_isi_violation_rate(spike_train_ms)
    ax.bar(x=xs, height=ys, width=bin_ms, color="gray", align="edge")
    ax.set_title(f'ISI violation rate ({isi_threshold_ms}ms): {rate*100:0.1f}%')
    ax.set_xlabel('time (ms)')

    ax = plt.subplot(n_rows, n_cols, 7)
    unit_extremum_waveforms = waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channels[unit_id]]
    ax.plot(sample_objects(unit_extremum_waveforms).T, label=unit_id, lw=0.5)
    ax.set_title(f'{unit_id} waveforms at ch {extremum_channels[unit_id]}')

    ax = plt.subplot(n_rows, n_cols, 8)
    unit_extremum_template = waveform_extractor.get_template(unit_id)[:, extremum_channels[unit_id]]
    ax.plot(unit_extremum_template.T, label=unit_id)
    ax.set_title(f'{unit_id} template at ch {extremum_channels[unit_id]}')

    ax = plt.subplot(n_rows, n_cols, 10)
    deoutlier_extremum_waveforms = remove_outlier(unit_extremum_waveforms)
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
        ax = plt.subplot(n_rows, 2, 9)
        unit_shank_waveform_by_session_adata = get_unit_shank_waveforms_by_session_adata(unit_id, channel_indices, sorting, waveform_extractor, extremum_channels, sessions, max_count_per_session=1000)
        for session_i in sessions.index:
            session_shank_waveform_adata = unit_shank_waveform_by_session_adata[unit_shank_waveform_by_session_adata.obs['session_i'] == session_i]
            ax.plot(session_shank_waveform_adata.X.mean(0), color=plt.cm.turbo(session_i / len(sessions)))
            
        ax = plt.subplot(1, n_cols, 3)
        unit_extremum_waveform_by_session_adata = get_unit_extremum_waveforms_by_session_adata(unit_id, sorting, waveform_extractor, extremum_channels, sessions, max_count_per_session=1000)
        trace_gap = -100
        for session_i in sessions.index:
            session_extremum_waveform_adata = unit_extremum_waveform_by_session_adata[unit_extremum_waveform_by_session_adata.obs['session_i'] == session_i]
            ax.plot(session_extremum_waveform_adata.X.mean(0) + session_i * trace_gap, color=plt.cm.turbo(session_i / len(sessions)))
        ax.set_yticks(np.arange(0, trace_gap * len(sessions), trace_gap), sessions['date'])
        ax.yaxis.tick_right()

        plt.savefig(f'{savepath}0.png', bbox_inches='tight') # Hack to remove 3D plot excessive margins.
        plt.close()

        plot_UMAP_by_session(sessions, unit_shank_waveform_by_session_adata, sessions.index, initdate, savepath=f'{savepath}1.png')

        left = plt.imread(f'{savepath}0.png')
        right = plt.imread(f'{savepath}1.png')
        right = right[:, right.shape[1]//3+100:-right.shape[1]//4-77]
        collate_lateral_figures(left, right)
        os.remove(f'{savepath}0.png')
        os.remove(f'{savepath}1.png')
        plt.savefig(f'{savepath}2.png', bbox_inches='tight')

        try:
            plot_isi_by_session(unit_id, sorting, sessions, sessions.index, initdate, savepath=f'{savepath}3.png')  
             
            left = plt.imread(f'{savepath}2.png')
            right = plt.imread(f'{savepath}3.png')
            collate_lateral_figures(left, right, split=(8, 5))
            os.remove(f'{savepath}2.png')
            os.remove(f'{savepath}3.png')

        except Exception as e:
            print(f'  Failed to save ISI\n{e}')     
            shutil.move(f'{savepath}2.png', savepath)    
            return 
    plt.savefig(savepath, bbox_inches='tight')
    plt.close()


def plot_traces(traces, sampling_frequency, channel_indices, title, savepath, trace_gap=250, shank_gap=500, fontsize=25):
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
    plt.ylabel(rf'{trace_gap} $\mu$V gap between traces')
    plt.savefig(savepath, bbox_inches='tight')
    plt.close()