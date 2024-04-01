import matplotlib 
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np

from probeinterface import generate_multi_columns_probe
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes

# region 1
A_active_channel_names = [
    f'A-0{channel_index:02d}' 
        for channel_index in [23, 7, 22, 6, 21, 5, 20, 4, 19, 3, 18, 2, 17, 1, 16, 0, 15, 31, 14, 30, 13, 29, 12, 28, 11, 27, 10, 26, 9, 25, 8, 24]
]
B1_active_channel_names = [
    f'B-0{channel_index:02d}' 
        for channel_index in [0, 15, 1, 14, 2, 13, 3, 12]
]

# region 2
B2_active_channel_names = [
    f'B-0{channel_idnex:02d}'
        for channel_idnex in [4, 11, 5, 10, 6, 9, 7, 8]
]
C_active_channel_names = [
    f'C-0{channel_index:02d}' 
        for channel_index in [17, 46, 19, 44, 21, 42, 23, 40, 25, 38, 27, 36, 29, 34, 31, 32, 33, 30, 35, 28, 37, 26, 39, 24, 41, 22, 43, 20, 45, 18, 47, 16]
]

active_channel_names = {
    'region1': A_active_channel_names + B1_active_channel_names, # e.g. CA1
    'region2': B2_active_channel_names + C_active_channel_names, # e.g. M1
}

channel_indices = np.array([
    [ 0,  1,  2,  3,  7, 6, 5, 4], #4,  5,  6,  7],
    [ 8,  9, 10, 11, 15, 14, 13, 12], #12, 13, 14, 15],
    [16, 17, 18, 19, 23, 22, 21, 20], #20, 21, 22, 23],
    [24, 25, 26, 27, 31, 30, 29, 28], #28, 29, 30, 31],
    [32, 33, 34, 35, 39, 38, 37, 36], #36, 37, 38, 39],
])
shank_locations = np.array([(0, 300), (100, 150), (200, 0), (300, 150), (400, 300)])
nrows, ncols = 4, 2
inter_electrode_distance = 30
electrode_radius = 10

def create_multi_shank_probe(savepath=None):
    plt.rcParams.update({'font.size': 8})
    n_shank = len(channel_indices)
    n_channel = channel_indices.size
    
    plt.figure(figsize=(20, 10))
    ax = plt.gca()

    probes = []
    for shank_channel_indices, shank_location in zip(channel_indices, shank_locations):
        probe = generate_multi_columns_probe(
            num_columns=ncols, num_contact_per_column=nrows, 
            xpitch=inter_electrode_distance, ypitch=inter_electrode_distance,
            contact_shapes='circle', contact_shape_params={'radius': electrode_radius}
        )
        probe.move(shank_location)
        probe.set_device_channel_indices(shank_channel_indices)

        plot_probe(probe, with_device_index=True, ax=ax)
        probes.append(probe)
        
    multi_shank_probe = combine_probes(probes)
    multi_shank_probe.set_device_channel_indices(channel_indices.flatten())

    plt.xlim(-150, 600)
    plt.ylim(-200, 600)
    plt.title(f'Probe - {n_channel}ch - {n_shank}shanks')
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    plt.close()
    return multi_shank_probe

def create_single_shank_probe(shank, savepath=None):
    plt.rcParams.update({'font.size': 8})
    
    plt.figure(figsize=(20, 10))
    ax = plt.gca()
    probe = generate_multi_columns_probe(
        num_columns=ncols, num_contact_per_column=nrows, 
        xpitch=inter_electrode_distance, ypitch=inter_electrode_distance,
        contact_shapes='circle', contact_shape_params={'radius': electrode_radius}
    )
    probe.set_device_channel_indices(np.argsort(channel_indices[shank]))
    plot_probe(probe, with_device_index=True, ax=ax)
    plt.xlim(-150, 150)
    plt.ylim(-200, 200)
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    plt.close()
    return probe