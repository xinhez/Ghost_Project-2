import spikeinterface.preprocessing as spre 

def preprocess(recording):
    recording = spre.bandpass_filter(recording, freq_min=300, freq_max=6000)    
    recording = spre.common_reference(recording, reference='global', operator='median')
    return recording