

plot_dir = '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/';
data_dir = '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/ExampleData/';


% load data 
load([data_dir 'cal_coeffs.mat']); % created by RAD_CAL_with_Linear_Fit_and_uncertainty_with_Monte_Carlo.m
load([data_dir 'pico_spectra.mat']); % QE Pro covering fluorescence range

IT = 160; % integration time of instrument in the field


% dark current correct
DN_E = spectra_struct.vectors(1, :) - spectra_struct.vectors(3, :);
DN_L = spectra_struct.vectors(2, :) - spectra_struct.vectors(4, :);


fh=figure;
hold
plot(spectra_struct.wvl, DN_L);
plot(spectra_struct.wvl, DN_E, 'r');
xlabel('Wavelength [nm]');
ylabel('DN')
title_str = 'DN corrected for DC';
title(title_str)
print_pdf(fh, plot_dir, title_str);


% % resample the example data to the instrument wavelengths we calibrated
% in the lab (we have just not got outside data for the system we got in
% the lab for this example)

DN_E_int = interp1(spectra_struct.wvl, DN_E, cal_coeffs.wvl,'linear','extrap');
DN_L_int = interp1(spectra_struct.wvl, DN_L, cal_coeffs.wvl,'linear','extrap');


% radiometric calibration

L = (DN_L_int/IT - cal_coeffs.offset)  ./ cal_coeffs.gain;




fh=figure;
hold
plot(cal_coeffs.wvl, L);
xlabel('Wavelength [nm]');
ylabel('L [W/m^2/sr/nm]')
title_str = 'Calibrated L';
title(title_str)
print_pdf(fh, plot_dir, title_str);

