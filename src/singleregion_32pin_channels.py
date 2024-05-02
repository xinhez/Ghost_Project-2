import matplotlib.pyplot as plt
import numpy as np

from probeinterface import generate_multi_columns_probe
from probeinterface.plotting import plot_probe
from probeinterface.utils import combine_probes


n_channel = 30

channel_indices = np.array([
    [ 5,  3,  1,  7,  9, 11], # 5
    [17, 15, 13, 19, 21, 23], # 4
    [29, 27, 25, 28, 26, 24], # 3
    [18, 20, 22, 16, 14, 12], # 2
    [ 6,  8, 10,  4,  2,  0], # 1
])
shank_locations = np.array([[0, 0], [150, 200], [300, 400], [450, 200], [600, 0]])
nrows, ncols = 3, 2
inter_electrode_distance = 50
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

    plt.xlim(-150, 800)
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
    probe.set_device_channel_indices(np.arange(len(channel_indices[shank])))
    plot_probe(probe, with_device_index=True, ax=ax)
    plt.xlim(-100, 150)
    plt.ylim(-200, 200)
    if savepath is not None:
        plt.savefig(savepath, bbox_inches='tight')
    plt.close()
    return probe