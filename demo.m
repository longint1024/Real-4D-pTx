clear
close all
clc
addpath('data/');
addpath('src/');

ID = [5]; %%% serial number of calibration

slab_thickness = 48; %%% defalt position: isocenter
%%% Users can manually move the geometric center of the calibration data to adjust the position of the slab
dwell_time = 10e-6;

%%% design Real-4D-pTx pulse
[rf,grad] = design_pTxSSWE(ID, slab_thickness, dwell_time);
%%% The calculation process is very slow. Please be patient and wait for around 20 minutes
%%% The default here is to balance passband ripple, stopband ripple, fat suppression, and RF energy, with pre-set hyperparameters. 
%%% If users are interested, they can manually adjust these hyperparameters, and there are no adjustment options provided in this demo.