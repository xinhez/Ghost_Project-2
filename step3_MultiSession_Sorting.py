import matplotlib
matplotlib.use('Agg')

import datetime
import glob
import matplotlib.pyplot as plt 
import numpy as np
import os
import pandas as pd
import spikeinterface.core as sc
import spikeinterface.curation as scu
import spikeinterface.extractors as se
import spikeinterface.widgets as sw
import sys 

from tqdm.auto import tqdm

sys.path.append('src')

from utils import surgery_dates, blackrock_channel_indices, plot_unit

if __name__ == '__main__':
    sorted_date = '20240112'
    mice = sorted([mouse_path.split(os.sep)[-1] for mouse_path in glob.glob(f'data{os.sep}sorted{os.sep}{sorted_date}{os.sep}**')])
    mice = ['6_2', '6_3', '6_7', '7_2']
    print(f'Plotting mice: {mice} sorted on {sorted_date}')


    for mouse in (pbar := tqdm(mice)):
        mouse_sorted_folder = f'data{os.sep}sorted{os.sep}{sorted_date}{os.sep}{mouse}'
        pbar.set_description(mouse)

        mouse_recording_si_path = f'{mouse_sorted_folder}{os.sep}processed'
        mouse_sorting_si_path = f'{mouse_sorted_folder}{os.sep}sorting'
        mouse_waveforms_si_path = f'{mouse_sorted_folder}{os.sep}waveforms'
        mouse_units_si_path = f'{mouse_sorted_folder}{os.sep}units'


        recording_processed = sc.load_extractor(mouse_recording_si_path)

        sorting = se.NpzSortingExtractor(f'{mouse_sorting_si_path}{os.sep}sorter_output{os.sep}firings.npz')
        sorting = scu.remove_excess_spikes(sorting, recording_processed) # spikeinterface https://github.com/SpikeInterface/spikeinterface/pull/1378

        waveform_extractor = sc.load_waveforms(
            folder=mouse_waveforms_si_path, with_recording=True, sorting=sorting
        )
        extremum_channels = sc.get_template_extremum_channel(waveform_extractor, peak_sign='neg')

        mouse_sessions = pd.read_csv(f'{mouse_sorted_folder}{os.sep}sessions.csv').sort_values(by='date')

        os.makedirs(mouse_units_si_path, exist_ok=True)
        for unit_id in sorting.unit_ids:
            pbar.set_description(f'{mouse} Plotting [unit {unit_id} / {len(sorting.unit_ids)}]')
            unit_plot_file = f'{mouse_units_si_path}{os.sep}{unit_id}.png'
            print(unit_plot_file)
            if True:# not os.path.isfile(unit_plot_file):
                plot_unit(waveform_extractor, extremum_channels, sorting, unit_id, blackrock_channel_indices, unit_plot_file, sessions=mouse_sessions)