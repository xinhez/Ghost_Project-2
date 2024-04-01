import anndata as ad 
import datetime
import matplotlib.pyplot as plt
import numpy as np
import scanpy
import spikeinterface.core as sc
import spikeinterface.widgets as sw
import sys 

sys.path.append('src')

from src.facts import *

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

def plot_unit(unit_id, session_info, init_date, waveform_extractors, channel_indices, savepath,
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
    n_rows = len(session_info['segment_path'].unique())
    n_cols = len(plot_types)
    plt.figure(figsize=(n_cols * subplot_size, n_rows * subplot_size))
    waveforms, templates, extremum_channels, dates = [], [], [], []
    for segment, segment_path in enumerate(session_info['segment_path'].unique()):
        waveform_extractor = waveform_extractors[segment]
        extremum_channel = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg')[unit_id]
        extremum_channel = np.where(waveform_extractors[segment].channel_ids == extremum_channel)[0].item()
        if len(channel_indices.shape) > 1:
            extremum_shank = get_shank(extremum_channel, channel_indices)
            extremum_channel_indices = channel_indices[extremum_shank]
        else:
            extremum_channel_indices = np.argsort(channel_indices)
        extremum_channels.append(extremum_channel)

        segment_waveforms = sample_objects(waveform_extractor.get_waveforms(unit_id)[:, :, extremum_channel_indices], max_n=100)
        segment_templates = waveform_extractor.get_template(unit_id)[:, extremum_channel_indices]
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

def plot_traces(traces, sampling_frequency, channel_indices, title, savepath, session_w=10, trace_gap=75, shank_gap=200, fontsize=25):
    n_shank, n_channel_per_shank = channel_indices.shape
    n_channel = channel_indices.size
    plt.rcParams.update({'font.size': fontsize})
    duration = traces.shape[1] / sampling_frequency / n_s_per_min
    plt.figure(figsize=(session_w * duration, 50))
    plt.title(f'{title} : {duration:0.2f} min')

    for shank_i, shank in enumerate(channel_indices):
        for channel_i, channel in enumerate(shank):
            y_baseline = trace_gap * (channel_i + shank_i * n_channel_per_shank) + shank_gap * shank_i
            
            plt.plot(traces[channel - channel_indices.min()]+y_baseline)
            plt.text(len(traces[channel - channel_indices.min()]), y_baseline - fontsize, f'ch{channel}')

    xticks_labels = list(range(round(duration) + 1))
    xticks_locs = [min * sampling_frequency * n_s_per_min for min in xticks_labels]
    plt.ylim(-trace_gap, n_channel * trace_gap + (n_shank - 1) * shank_gap)
    plt.xticks(ticks=xticks_locs, labels=xticks_labels)
    plt.xlabel('min')
    plt.ylabel(rf'{trace_gap} $\mu$V gap between traces')
    plt.savefig(savepath, bbox_inches='tight')
    plt.close()
    
