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
    'M9_6':  '20240128',
    'M10_8': '20240128',
    'M11_4': '20240128',
    'M9_4':  '20240211',
    'M9_7':  '20240211',
    'M10_1': '20240211',
    'M10_5': '20240211',
    'M10_6': '20240222',
    'M9_8':  '20240316',
    'M15_2': '20240316',
    'M15_5': '20240316',
    'M16_1': '20240316',
    'M16_2': '20240316',
    'M16_6': '20240316',
    'M17_2': '20240316',
    'M15_3': '20240317',
    'M15_7': '20240317',
    'M16_7': '20240317',
    'M17_5': '20240317',
    # Deep Brain 
    'D12_6': '20240304',
    'D13_4': '20240304',
    'D13_8': '20240304',
    'D14_4': '20240304',
    'D14_6': '20240304',
}

n_day_per_week = 7
n_s_per_min = 60
n_ms_per_s = 1000
ms_before, ms_after = 2, 2
total_week = 12

window_ms = 150
bin_ms = 1.5
isi_threshold_ms = 1.5
detect_interval_ms = 0.33