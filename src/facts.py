probe_designs = {
    # Multiple Region
    'M9_7': 'multiregion_80pin',
    'M10_6': 'multiregion_80pin',
    # Lateral Implantation
    'L11_9': 'multiregion_80pin',
    'L16_8': 'multiregion_80pin',
    'L17_7': 'multiregion_80pin',
    'L14_5': 'singleregion_80pin',
    # Long Term
    '1_5': 'singleregion_32pin',
    '5_7': 'singleregion_32pin',
    '6_2': 'singleregion_32pin',
    '6_7': 'singleregion_32pin',
    '7_2': 'singleregion_32pin',
}

surgery_dates = {
    # Lateral Implantation
    'L11_9': '20240320',
    'L14_5': '20240320',
    'L16_8': '20240320',
    'L17_7': '20240320',
    # Long Term
    '1_5': '20230627',
    '5_7': '20230805',
    '6_2': '20231019',
    '6_7': '20231016',
    '7_2': '20231019',
}

n_day_per_week = 7
n_s_per_min = 60
n_ms_per_s = 1000
ms_before, ms_after = 2, 2
total_week = 20

window_ms = 150
bin_ms = 1.5
isi_threshold_ms = 1.5