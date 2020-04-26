# LFP-Type2-and-Type3-Detection
Detecting Type 2 and 3 events from LFP data: This is code for the detection of Type 2 and 3 events for the manuscript titled "Microscale physiological events on the human cortical surface" by Paulk et al., preprint located at https://www.biorxiv.org/content/10.1101/770743v1  

The criteria for detecting Type 2 and 3 events are as follows:
1) absolute waveforms were detected which were >25 µV in amplitude with a series of steps to clean the waveforms and remove noise
2) detected waveforms had a correlation above 0.8 with a template event waveform
3) the second derivative at the onset of the recording was greater than 2
4) the voltage in the 100 ms preceding the event onset is less than 25 µV to reduce the chances of capturing oscillations


Input data:
LFP data should be in the form of channel x sample, sampled at 1000 Hz, low pass filtered at 500 Hz, with the voltage scale in microVolts.
