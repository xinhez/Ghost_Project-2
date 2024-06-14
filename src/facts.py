probe_designs = {
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
    # Multiple Region
    'M9_4': 'multiregion_80pin',
    'M9_6': 'multiregion_80pin',
    'M9_7': 'multiregion_80pin',
    'M9_8': 'multiregion_80pin',
    'M10_1': 'multiregion_80pin',
    'M10_5': 'multiregion_80pin',
    'M10_6': 'multiregion_80pin',
    'M10_8': 'multiregion_80pin',
    'M11_4': 'multiregion_80pin',
    'M15_2': 'multiregion_80pin',
    'M15_3': 'multiregion_80pin',
    'M15_5': 'multiregion_80pin',
    'M15_7': 'multiregion_80pin',
    'M16_1': 'multiregion_80pin',
    'M16_2': 'multiregion_80pin',
    'M16_6': 'multiregion_80pin',
    'M16_7': 'multiregion_80pin',
    'M17_2': 'multiregion_80pin',
    'M17_5': 'multiregion_80pin',
}

surgery_dates = {
    # Lateral Implantation
    'L11_9': '20240320',
    'L14_5': '20240320',
    'L16_8': '20240320',
    'L17_7': '20240320',
    # BlackRock Long Term
    '1_5': '20230627',
    # Multi Region
    'M9_6': '240128',
    'M10_8': '240128',
    'M11_4': '240128',
    'M9_4': '240211',
    'M9_7': '240211',
    'M10_1': '240211',
    'M10_5': '240211',
    'M10_6': '240222',
    'M9_8': '240316',
    'M15_2': '240316',
    'M15_5': '240316',
    'M16_1': '240316',
    'M16_2': '240316',
    'M16_6': '240316',
    'M17_2': '240316',
    'M15_3': '240317',
    'M15_7': '240317',
    'M16_7': '240317',
    'M17_5': '240317',
    # Deep Brain 
    'D12_6': '240304',
    'D13_4': '240304',
    'D13_8': '240304',
    'D14_4': '240304',
    'D14_6': '240304',
}

n_day_per_week = 7
n_s_per_min = 60
n_ms_per_s = 1000
ms_before, ms_after = 2, 2
total_week = 20

window_ms = 150
bin_ms = 1.5
isi_threshold_ms = 1.5