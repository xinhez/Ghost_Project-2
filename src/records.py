import numpy as np

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

longterm_segments = {
    '1_5': {
        '20230707': [(0, 20)],
        '20230714': [(0, 20)],
        '20230721': [(0, 20)],
        '20230727': [(0, 20)],
        '20230821': [(0, 20)],
        '20230831': [(10, 30)],
        '20230921': [(10, 30)],
        '20231006': [(0, 20)],
        '20231011': [(0, 20)],
        '20231018': [(0, 15), (25, 30)],
        '20231027': [(10, 30)],
        '20231102': [(0, 20)],
        '20231122': [(0, 20)],
        '20231207': [(0, 20)],
        '20231214': [(5, 20), (25, 30)],
        '20231228': [(6, 15), (16, 18), (21, 30)],
        '20240103': [(0, 20)],
        '20240117': [(10, 30)],
    }
}

blackrock_channel_indices = np.array([ 
    [ 5,  3,  1,  7,  9, 11], 
    [17, 15, 13, 19, 21, 23], 
    [29, 27, 25, 28, 26, 24], 
    [18, 20, 22, 16, 14, 12], 
    [ 6,  8, 10,  4,  2,  0],
])

intan_channel_indices = np.array([
    # [21, 22, 23, 20, 19, 18],# [24, 23, 22, 27, 26, 25], 
    # [15, 16, 17, 14, 13, 12],# [ 0, 29, 28,  3,  2,  1], 
    # [ 9, 10, 11,  6,  5,  4],# [ 6,  5,  4,  9,  8,  7], 
    # [ 1,  2,  3,  0, 31, 30],# [12, 11, 10, 15, 14, 13], 
    # [27, 28, 29, 26, 25, 24],# [18, 17, 16, 21, 20, 19],
    [19, 20, 21, 18, 17, 16],
    [13, 14, 15, 12, 11, 10],
    [ 7,  8,  9,  6,  5,  4],
    [ 1,  2,  3,  0, 29, 28],
    [25, 26, 27, 24, 23, 22]
])

shank_locations = np.array([[0, 0], [150, 200], [300, 400], [450, 200], [600, 0]])

n_day_per_month= 30
n_s_per_min = 60
n_ms_per_s = 1000

window_ms, bin_ms = 100, 1.5
isi_threshold_ms = 1.5

memory_limit = '20G'
min_recording_duration = 10

sorter_parameters = {
    'detect_sign': -1,
    'adjacency_radius': 120, 
    'freq_min': 300, 
    'freq_max': 3000,
    'filter': True,
    'whiten': True,  
    'clip_size': 50,
    'detect_threshold': 4,
    'detect_interval': 3, # 0.3ms 
    'tempdir': 'sorter_tmp',
}

sorter5_parameters = {
    "scheme": "2",  # '1', '2', '3'
    # "scheme": "Which sorting scheme to use: '1, '2', or '3'",
    "detect_threshold": 5.5,  # this is the recommended detection threshold
    # "detect_threshold": "Detection threshold - recommend to use the default",
    "detect_sign": -1,
    # "detect_sign": "Use -1 for detecting negative peaks, 1 for positive, 0 for both",
    "detect_time_radius_msec": 0.5,
    # "detect_time_radius_msec": "Determines the minimum allowable time interval between detected spikes in the same spatial region",
    "snippet_T1": 20,
    # "snippet_T1": "Number of samples before the peak to include in the snippet",
    "snippet_T2": 20,
    # "snippet_T2": "Number of samples after the peak to include in the snippet",
    "npca_per_channel": 3,
    # "npca_per_channel": "Number of PCA features per channel in the initial dimension reduction step",
    "npca_per_subdivision": 10,
    # "npca_per_subdivision": "Number of PCA features to compute at each stage of clustering in the isosplit6 subdivision method",
    "snippet_mask_radius": 250,
    # "snippet_mask_radius": "Radius of the mask to apply to the extracted snippets",
    "scheme1_detect_channel_radius": -1,
    # "scheme1_detect_channel_radius": "Channel radius for excluding events that are too close in time in scheme 1",
    "scheme2_phase1_detect_channel_radius": 200,
    # "scheme2_phase1_detect_channel_radius": "Channel radius for excluding events that are too close in time during phase 1 of scheme 2",
    "scheme2_detect_channel_radius": 50,
    # "scheme2_detect_channel_radius": "Channel radius for excluding events that are too close in time during phase 2 of scheme 2",
    "scheme2_max_num_snippets_per_training_batch": 200,
    # "scheme2_max_num_snippets_per_training_batch": "Maximum number of snippets to use in each batch for training during phase 2 of scheme 2",
    "scheme2_training_duration_sec": 60 * 5,
    # "scheme2_training_duration_sec": "Duration of training data to use in scheme 2",
    "scheme2_training_recording_sampling_mode": "uniform",
    # "scheme2_training_recording_sampling_mode": "initial or uniform",
    "scheme3_block_duration_sec": 60 * 30,
    # "scheme3_block_duration_sec": "Duration of each block in scheme 3",
    "freq_min": 300,
    # "freq_min": "High-pass filter cutoff frequency",
    "freq_max": 3000,
    # "freq_max": "Low-pass filter cutoff frequency",
    "filter": True,
    # "filter": "Enable or disable filter",
    "whiten": True,  # Important to do whitening
    # "whiten": "Enable or disable whitening",
    "temporary_base_dir": 'sorter_temp',
    # "temporary_base_dir": "Temporary directory base directory for storing cached recording",
    "n_jobs_for_preprocessing": -1,
    # "n_jobs_for_preprocessing": "Number of parallel jobs for creating the cached recording",
}

ms_before, ms_after = 1, 2