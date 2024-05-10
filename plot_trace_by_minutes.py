import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import spikeinterface.extractors as se 
import spikeinterface.preprocessing as spre
import sys

from tqdm.auto import tqdm
sys.path.append('src')

from src.importrhdutilities import load_file

subjects = ['D13_4', 'D14_6', 'D12_6', 'D13_8']

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

sortdate = '240507'
n_s_per_min = 60

def plot_traces(traces, sampling_frequency, channel_indices, title, savepath, session_w=10, trace_gap=75, shank_gap=150, fontsize=25):
    n_shank, n_channel_per_shank = channel_indices.shape
    n_channel = channel_indices.size
    plt.rcParams.update({'font.size': fontsize})
    duration = traces.shape[1] / sampling_frequency / n_s_per_min
    plt.figure(figsize=(session_w * duration, 50))
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

for subject in (pbar := tqdm(subjects)):
    session_path = glob.glob(f'data/raw/{sortdate}*/{subject}*')[0]
    
    recording_paths = sorted(glob.glob(f'{session_path}/*.rhd'))

    for recording_index, recording_path in enumerate(recording_paths):
        pbar.set_description(f'{subject} : {recording_index + 1}/{len(recording_paths)}')

        trace_file = recording_path.replace('.rhd', '.png').replace('raw', 'processed')
        if os.path.isfile(trace_file): continue

        output_folder = '/'.join(trace_file.split('/')[:-1])
        
        os.makedirs(output_folder, exist_ok=True)

        raw_data, data_present = load_file(recording_path)
        if data_present:
            recording_channel_names = [channel_info['native_channel_name'] for channel_info in raw_data['amplifier_channels']]
            sampling_frequency = raw_data['frequency_parameters']['amplifier_sample_rate']
            active_channel_indices = [recording_channel_names.index(active_channel_name) for active_channel_name in active_channel_names]
            traces = raw_data['amplifier_data'][active_channel_indices]

            recording = se.NumpyRecording(traces_list=traces.T, sampling_frequency=sampling_frequency)
            recording = spre.bandpass_filter(recording, freq_min=300, freq_max=3000)
            recording = spre.common_reference(recording, reference='global', operator='median')

            traces = recording.get_traces().T
            plot_traces(traces, sampling_frequency, channel_indices, recording_path.split('/')[-1], trace_file)