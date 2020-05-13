# LFP-Type2-and-Type3-Detection
Detecting Type 2 and 3 events from LFP data: This is code for the detection of Type 2 and 3 events for the manuscript titled "Microscale physiological events on the human cortical surface" by Paulk et al., preprint located at https://www.biorxiv.org/content/10.1101/770743v1  

The criteria for detecting Type 2 and 3 events are as follows:
1) absolute waveforms were detected which were >25 µV in amplitude with a series of steps to clean the waveforms and remove noise
2) detected waveforms had a correlation above 0.8 with a template event waveform
3) the second derivative at the onset of the recording was greater than 2
4) the voltage in the 100 ms preceding the event onset is less than 25 µV to reduce the chances of capturing oscillations

The example .m file includes comments and step-by-step process of detecting the waveforms, cleaning the waveforms, correlating them to the templates, and finally plotting the waveforms.  

The MATLAB code is currently written for a single channel of LFP input. It was tested in MATLAB 2019a and 2016b.

Example data is included (ExampleData.zip) in a .mat file (Examp.mat), with the variables datalfp and sampling rate (DecFS)

Since there can be too-large deflections in the data set, some of the steps involve removing large amplitude deflections or sharp transitions.

Input data:
LFP data should be in the form of channel x sample, sampled at 1000 Hz, low pass filtered at 500 Hz, with the voltage scale in microVolts.

Output: 
Figure of the snippet waveforms as well as a stem raster plot of the waveform times in the recording.
