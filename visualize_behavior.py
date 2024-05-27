import matplotlib.pyplot as plt
import numpy as np 
import os
import pandas as pd
import seaborn as sns
import spikeinterface.core as sc
import spikeinterface.curation as scu
import spikeinterface.extractors as se

from tqdm.auto import tqdm 

s_before_cue = 1 
s_after_cue = 2
n_segment = 14
n_ms_per_s = 1000
bin_size_ms = 100
n_total_trial = 100

for threshold in [3.5, 4.5, 5.5, 6.5]:
    for subject in (pbar := tqdm(['M15_2', 'M15_3', 'M15_5', 'M15_7', 'M16_1'])):
        session_info = pd.read_csv(f'data/processed/{subject}/info.csv')
        n_segment = len(session_info)
        for region in ['CA1', 'M1']:
            pbar.set_description(f'threshold {threshold} {subject} - {region}')
            recordings = [sc.load_extractor(f'data/processed/{subject}/{region}/recordings/segment{segment_index}') for segment_index in range(n_segment)]
            recording = sc.concatenate_recordings(recordings)

            sorting = se.NpzSortingExtractor(f'data/processed/{subject}/{region}/sortings-{threshold}/sorter_output/firings.npz')
            sorting = scu.remove_excess_spikes(sorting, recording)
            sortings = sc.split_sorting(sorting, recordings)
            sortings = [sc.select_segment_sorting(sortings, segment_indices=segment) for segment in range(len(recordings))]

            heatmap_folder = f'data/processed/{subject}/{region}/heatmaps-{threshold}'
            os.makedirs(heatmap_folder, exist_ok=True)

            n_frames_per_ms = int(sorting.sampling_frequency / n_ms_per_s)
            for unit_id in sorting.unit_ids:
                pbar.set_description(f'threshold {threshold} {subject} - {region} - {unit_id} / {len(sorting.unit_ids)}')

                if os.path.isfile(f'{heatmap_folder}/{unit_id}.png'):
                    continue
                trial_bins = []
                for segment_index in range(n_segment):
                    segment_info = session_info[session_info['segment_index']==segment_index]

                    controller_date = segment_info['segment_path'].item().split('_')[-2]
                    controller_date = controller_date[2:] + controller_date[:2]
                    controller_file = f'data/raw/controller/{controller_date}/{subject.replace("_", "-")}.TXT'
                    events = pd.read_csv(controller_file)
                    # If there are multiple headers that day, the last header marks the start of the actual experiment. This might happen when we use the subject RFID tag duing the equipment check.
                    state_ids = np.where(events['state'] == 'state')[0]
                    if len(state_ids) > 0:
                        events = events.iloc[state_ids[-1]+1:].reset_index(drop=True)
                    s_trials, f_trials = [], []
                    for trial in range(n_total_trial):
                        trial_indices = np.where(events['trial'] == trial)[0]
                        if any(events['state'][trial_indices] == 1):
                            s_trials.append(trial)
                        else:
                            f_trials.append(trial)
                    assert len(s_trials) + len(f_trials) == n_total_trial

                    trial_starts = eval(segment_info['trial_starts'].item())
                    trial_ends = eval(segment_info['trial_ends'].item())

                    spike_train = sortings[segment_index].get_unit_spike_train(unit_id)

                    segment_trial_bins = []
                    for trial, (trial_start, trial_end) in enumerate(zip(trial_starts, trial_ends)):
                        if trial not in s_trials: continue

                        t_start = trial_start - s_before_cue * sorting.sampling_frequency
                        t_end = trial_start + s_after_cue * sorting.sampling_frequency
                        trial_spikes = spike_train[(spike_train >= t_start) & (spike_train < t_end)] - trial_start
                        
                        segment_trial_bins.append(np.histogram(trial_spikes, bins=np.arange(-s_before_cue*sorting.sampling_frequency, s_after_cue*sorting.sampling_frequency+1, n_frames_per_ms*bin_size_ms))[0])
                    segment_trial_bins = np.array(segment_trial_bins)
                    trial_bins.append(segment_trial_bins.mean(0))
                trial_bins = np.array(trial_bins)
                plt.figure()
                sns.heatmap(trial_bins)
                plt.title(f'{subject} - {region} - {unit_id}')
                plt.yticks(np.arange(0, n_segment)+0.5, np.arange(1, n_segment+1))
                plt.ylabel('session')
                plt.xticks(np.arange(30), [f'{number:0.1f}' for number in np.arange(-s_before_cue, s_after_cue, 0.1)], rotation=90)
                plt.xlabel('cue onsite (s)')
                plt.savefig(f'{heatmap_folder}/{unit_id}.png')
                plt.close()