%% Radiometric Calibration using Linear Fit and 
%% using Monte Carlo for Propagation of L and DN Uncertainty 
%% to Gains and Offsets
%
% ahueni, 2019
% 
% Prepared for the SENSECO COST action
%
%
%%

plot_dir = '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/';
data_dir = '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/ExampleData/';

Linewidth = 3; % used for plotting

%% Get Data

load([data_dir 'DN_L_CAL.mat'])
load([data_dir 'L_Sphere.mat'])
load([data_dir 'uL.mat'])
load([data_dir 'STD_DN.mat'])


wvl_sphere = L_Sphere(:,1);
L_sphere = L_Sphere(:,2:end);

u_sphere = u_L(:,2); % relative uncertainty at k = 2
u_sphere = u_sphere / 2; % set to k=1
u_sphere_wvl = u_L(:,1); 

u_sphere_long_term_temporal_stability = 0.4; % relative uncertainty

IT = 15; % integration time in milliseconds

wvl_instr = DN_L_CAL(:,1); % wavelengths of instrument

DN_l_dc = DN_L_CAL(:,2:end); % mean values per intensity corrected for dark current

DN_l_dc = DN_l_dc / IT; % normalise by integration time 

% absolute uncertainty, taken from measurement noise divided by square root of N = 100 -> 10

u_DN = zeros(length(wvl_sphere), 4);
u_DN(:,1) = StdDeviation_5 ./ 10;
u_DN(:,2) = StdDeviation_100 ./ 10;
u_DN(:,3) = StdDeviation_1000 ./ 10;
u_DN(:,4) = StdDeviation_10000 ./ 10;

u_DN = u_DN / IT;  % normalise by integration time 


fh=figure;
plot(wvl_sphere, u_DN)
legend('fL5', 'fL100', 'fL1000', 'fL10000')
xlabel('Wavelength [nm]');
ylabel('u(DN)')
title_str = 'DN Uncertainties estimated from standard deviation';
title(title_str)
print_pdf(fh, plot_dir, title_str);


% plot input data
fh=figure;
semilogy(wvl_sphere, L_sphere);
xlabel('Wavelength [nm]');
ylabel('log(L [W/m^2/sr/nm])')
title_str = 'Sphere Radiances';
title(title_str)
print_pdf(fh, plot_dir, title_str);

fh=figure;
plot(wvl_instr, DN_l_dc);
xlabel('Wavelength [nm]');
ylabel('DN DC_Corr')
title_str = 'DNs (DC corrected)';
title(title_str)
print_pdf(fh, plot_dir, title_str);


fh=figure;
plot(wvl_sphere, L_sphere(:,1));
xlabel('Wavelength [nm]');
ylabel('L [W/m^2/sr/nm]')
title_str = 'Sphere Radiance: lowest intensity';
title(title_str)
print_pdf(fh, plot_dir, title_str);



% adjust lamp wavelengths to instrument wavelengths by resampling (ideally, do spectral convolution here)

L_sphere_interp = zeros(size(DN_l_dc));

for i =1: size(L_sphere, 2)
    L_sphere_interp(:,i) = interp1(wvl_sphere, L_sphere(:,i), wvl_instr);
end

% adjust sphere uncertainty wavelengths to instrument wavelengths by
% resampling of u_sphere
u_sphere = interp1(u_sphere_wvl, u_sphere, wvl_instr);


%% Get gain and offset (instrument sensitivity) using linear interpolation: this is the regular radiometric calibration

band_to_show = 20;

gain = zeros(size(wvl_instr));
offset = zeros(size(wvl_instr));

% loop over all bands and do a first order fit
for i=1:size(DN_l_dc,1)

    c = polyfit(L_sphere(i,:), DN_l_dc(i,:), 1);
    
    gain(i) = c(1);
    offset(i) = c(2);

    if i == band_to_show
        % show a plot of the fit for a selected band
        
        fh=figure;
        hold
        plot(L_sphere(i,:), DN_l_dc(i,:), 'o', 'MarkerSize', 12, 'MarkerFace', 'b')
        
        x = linspace(min(L_sphere(i,:)),max(L_sphere(i,:)));
        y=polyval(c,x);
        
        plot(x, y, 'r', 'Linewidth', Linewidth)
        
        xlabel('Radiance');
        ylabel('mean(DN dc_corr)')
        title_str = 'Example of linear fit to get gain and offset';
        title(title_str)
        
        text(round(mean(x)), round(mean(y)*0.9), ['gain =' num2str(c(1)) ', offset = ' num2str(c(2))])
        
        print_pdf(fh, plot_dir, title_str);

    end

end

cal_coeffs.gain = gain;
cal_coeffs.offset = offset;
cal_coeffs.wvl = wvl_instr;
save([data_dir 'cal_coeffs'], 'cal_coeffs');



fh=figure;
plot(wvl_instr, gain);
xlabel('Wavelength [nm]');
ylabel('Gain (DN/(W/m^2/sr/nm)')
title_str = 'Responsivity (Gain) of instrument';
title(title_str)
print_pdf(fh, plot_dir, title_str);

fh=figure;
plot(wvl_instr, offset);
xlabel('Wavelength [nm]');
ylabel('Offset (DN)')
title_str = 'Offset - of instrument';
title(title_str)
print_pdf(fh, plot_dir, title_str);


%% Get Calibration Coefficients using linear interpolation

band_to_show = 20;

gain_coeff = zeros(size(wvl_instr));
offset_coeff = zeros(size(wvl_instr));

% loop over all bands and do a first order fit
for i=1:size(DN_l_dc,1)

    c = polyfit(DN_l_dc(i,:), L_sphere(i,:), 1);
    
    gain_coeff(i) = c(1);
    offset_coeff(i) = c(2);

    if i == band_to_show
                
        fh=figure;
        hold
        plot(DN_l_dc(i,:), L_sphere(i,:), 'o', 'MarkerSize', 12, 'MarkerFace', 'b')
        
        x = linspace(min(DN_l_dc(i,:)),max(DN_l_dc(i,:)));
        y=polyval(c,x);
        
        plot(x, y, 'r', 'Linewidth', Linewidth)
        
        xlabel('mean(DN dc-corr)');
        ylabel('Radiance')
        title_str = 'Example of linear fit to get rad cal coefficients';
        title(title_str)
        
        text(round(mean(x)), mean(y)*0.8, ['c1 =' num2str(c(1)) ', c0 = ' num2str(c(2))])
        
        print_pdf(fh, plot_dir, title_str);
        
        
        % zoomed version on lower 3 data points
        xlim([min(DN_l_dc(i,:)), max(DN_l_dc(band_to_show,3))*1.1])
        title_str = 'Example of linear fit to get rad cal coefficients - zoom';
        title(title_str)
        
        print_pdf(fh, plot_dir, title_str);

    end

end




fh=figure;
plot(wvl_instr, gain_coeff);
xlabel('Wavelength [nm]');
ylabel('Rad Coeff (c1) (DN/(W/m^2/sr/nm)')
title_str = 'Rad Coeff (c1) of instrument';
title(title_str)
print_pdf(fh, plot_dir, title_str);

fh=figure;
plot(wvl_instr, offset_coeff);
xlabel('Wavelength [nm]');
ylabel('c0 (DN)')
title_str = 'Rad Coeff (c0) of instrument';
title(title_str)
print_pdf(fh, plot_dir, title_str);



%% Data radiometric calibration example

for i = 1:size(DN_l_dc, 2)

    L_calibrated = (DN_l_dc(:,i) - offset)  ./ gain;

    fh=figure;
    plot(wvl_instr, L_calibrated);
    hold
    plot(wvl_instr, L_sphere(:,i), 'r');
    legend('L Calibrated', 'L Sphere');
    xlabel('Wavelength [nm]');
    ylabel('L [W/m^2/sr/nm]')
    title_str = ['Calibrated Sphere Radiance - Intensity ' num2str(i)];
    title(title_str)
    print_pdf(fh, plot_dir, title_str);

end


%% combine uncertainties of sphere as relative uncertainties: 
u_sphere_combo = (u_sphere.^2 + u_sphere_long_term_temporal_stability.^2).^(0.5);


fh=figure;
plot(wvl_instr, u_sphere, 'r')
hold
plot(wvl_instr, ones(size(wvl_instr))*u_sphere_long_term_temporal_stability, 'g')
plot(wvl_instr, u_sphere_combo, 'm')

legend(escape_('u_sphere'), escape_('u_sphere_long_term_temporal_stability'), escape_('u_sphere_combo'))

yl = ylim;
ylim([0, yl(2)])
xlabel('Wavelength [nm]');
ylabel('relative uncertainty [%]')
title_str = ['Sources of Sphere Uncertainty and Combined Uncertainty'];
title(title_str)
print_pdf(fh, plot_dir, title_str);

%% ---------------------------------
%% MONTE CARLO follows hereafter ...
%% ---------------------------------

%% prepare realisations

N = 100;

L_realisations = zeros(size(L_sphere, 2), size(L_sphere, 1), N);

% Create Radiances realisations: loop over all intensities
for i=1:size(L_sphere, 2)
    
    random_numbers = randn(N, 1); % same random number for all spectral bands
    
    u_sphere_combo_abs = u_sphere_combo / 100 .* L_sphere(:,i);
    
    [L_realisations(i,:,:), ~] = get_realisations_gauss_dist(N, L_sphere(:,i), u_sphere_combo_abs, random_numbers);
    
    % plot the realisations per intensity
    fh=figure;
    plot(wvl_instr, squeeze(L_realisations(i,:,:)))
    
    xlabel('Wavelength [nm]');
    ylabel('L [W/m^2/sr/nm]')
    title_str = ['Sphere Radiance Realisations - Intensity ' num2str(i)];
    title(title_str)
    print_jpeg(fh, plot_dir, title_str);
    
    % plot an example of the sphere radiance, a selected realisation and
    % their error vector
    fh=figure;
    subplot(2,1,1);
    hold
    plot(wvl_sphere, L_sphere(:,i));
    plot(wvl_sphere, L_realisations(i,:,1), 'r');
    xlabel('Wavelength [nm]');
    ylabel('L [W/m^2/sr/nm]')
    title_str = 'Sphere Radiance: lowest intensity and example realisation ';
    legend('mean(L)', 'L realisation')
    title(title_str)
    subplot(2,1,2);
    error_vector = squeeze(L_realisations(i,:,1))' - L_sphere(:,i);
    plot(wvl_sphere, error_vector);
    xlabel('Wavelength [nm]');
    ylabel('L [W/m^2/sr/nm]')
    title_str = 'Sphere Radiance Simulated Error';
    title(title_str)    
    title_str = ['L Sphere: intensity ' num2str(i) ' and example realisation '];
    print_pdf(fh, plot_dir, title_str);

    
    

end


% plot showing all realisations

fh=figure;
plot(wvl_sphere, squeeze(L_realisations(:,:,1)), 'Linewidth', Linewidth)
legend('fL5', 'fL100', 'fL1000', 'fL10000')
xlabel('Wavelength [nm]');
ylabel('L [W/m^2/sr/nm]')
title_str = 'Sphere Radiance Realisation Example for all intensities';
title(title_str)
print_jpeg(fh, plot_dir, title_str);







% Create DN realisations

DN_realisations = zeros(size(L_sphere, 2), size(L_sphere, 1), N);

for i=1:size(L_sphere, 2)
        
    [DN_realisations(i,:,:), ~] = get_realisations_gauss_dist(N, DN_l_dc(:,i), u_DN(:,i));
    
    fh=figure;
    plot(wvl_instr, squeeze(DN_realisations(i,:,:)))
    
    xlabel('Wavelength [nm]');
    ylabel('mean(DN dc-corr)')
    title_str = ['DN Realisations - Intensity ' num2str(i)];
    title(title_str)
    print_jpeg(fh, plot_dir, title_str);
    
    fh=figure;
    subplot(2,1,1);
    hold
    plot(wvl_sphere, DN_l_dc(:,i));
    plot(wvl_sphere, DN_realisations(i,:,1), 'r');
    xlabel('Wavelength [nm]');
    ylabel('DN')
    title_str = ['DNs: intensity ' num2str(i) ' and example realisation '];
    legend('mean(DN)', 'DN realisation')
    title(title_str)
    subplot(2,1,2);
    error_vector = squeeze(DN_realisations(i,:,1))' - DN_l_dc(:,i);
    plot(wvl_sphere, error_vector);
    xlabel('Wavelength [nm]');
    ylabel('DN')
    title_str = 'DN Simulated Error';
    title(title_str)    
    title_str = ['DNs: intensity ' num2str(i) ' and example realisation '];
    print_pdf(fh, plot_dir, title_str);
    
    

end


Linewidth = 3;

fh=figure;
plot(wvl_sphere, squeeze(DN_realisations(:,:,1)), 'Linewidth', Linewidth)
legend('fL5', 'fL100', 'fL1000', 'fL10000')
xlabel('Wavelength [nm]');
ylabel('DN')
title_str = 'DN Realisation Example for all intensities';
title(title_str)
print_jpeg(fh, plot_dir, title_str);



%% Monte Carlo Simulation
% Only run simulation if required; otherwise load saved simulation results

run_sim = false

if run_sim
    
    gains = zeros(N^2, length(wvl_instr));
    offsets = zeros(N^2, length(wvl_instr));
    cnt = 1;
    
    progressbar_h = progressbar( [],0,['Calculating realisations for gain/offset fits' '...' ]);
    
    for n=1:N
               
        for m = 1:N
                       
            for i=1:size(DN_l_dc,1)
                
                c = polyfit(L_realisations(:,i,n), DN_realisations(:,i,m), 1);               
                gains(cnt,i) = c(1);
                offsets(cnt,i) = c(2);

            end
            
            progressbar( progressbar_h,1/(N^2));
            
            cnt = cnt + 1;
            
        end        
        
    end
    
    progressbar( progressbar_h,-1 );
    
    save([data_dir 'gain_realisations.mat'], 'gains');
    save([data_dir 'offset_realisations.mat'], 'offsets');
    
else
    
    load([data_dir 'gain_realisations.mat']);
    load([data_dir 'offset_realisations.mat']);
       
end


fh=figure;
plot(wvl_instr, gains);
xlabel('Wavelength [nm]');
ylabel('Gain (DN/(W/m^2/sr/nm)')
title_str = 'Responsivity (Gain) Realisations of instrument';
title(title_str)
print_jpeg(fh, plot_dir, title_str);

fh=figure;
plot(wvl_instr, offsets);
xlabel('Wavelength [nm]');
ylabel('Offset (DN)')
title_str = 'Offset Realisations of instrument';
title(title_str)
print_jpeg(fh, plot_dir, title_str);





if run_sim
    
    % get uncertainty of gain and offset

    u_gain = std(gains);
    u_offset = std(offsets);

    % calculat correlation between gain and offset and their combined uncertainty
    R = corrcoef(gain, offset);
    u_gain_offset = ((u_gain) .* (u_offset) .* R(2,1))';


    % save uncertainties
    u_rad_coeffs.u_gain = u_gain;
    u_rad_coeffs.u_offset = u_offset;
    u_rad_coeffs.u_gain_offset = u_gain_offset;    
    
    save([data_dir 'u_rad_coeffs.mat'], 'u_rad_coeffs');
else
    load([data_dir 'u_rad_coeffs.mat']);
end



fh=figure;
plot(wvl_instr, u_rad_coeffs.u_gain);
xlabel('Wavelength [nm]');
ylabel('u(Gain) (DN/(W/m^2/sr/nm)')
title_str = 'u(Gain) of instrument';
title(title_str)
print_pdf(fh, plot_dir, title_str);

fh=figure;
plot(wvl_instr, u_rad_coeffs.u_offset);
xlabel('Wavelength [nm]');
ylabel('Offset (DN)')
title_str = 'u(Offset) of instrument';
title(title_str)
print_pdf(fh, plot_dir, title_str);

fh=figure;
plot(gain, offset, 'o')
xlabel('gain');
ylabel('offset');
title_str = ['Correlation between gain and offset: R=' num2str(R(2,1))];
title(title_str)
print_pdf(fh, plot_dir, title_str);


%% Application to measured L


i = size(DN_l_dc, 2); % select the highest intensity as show-case

L_calibrated = (DN_l_dc(:,i) - offset)  ./ gain; % apply calibration coefficients to transform DN to radiance

uDN_ = u_DN(:,i); % the uncertainty of the DN's used for this calculation, taken from the uncertainties calculated from the DNs measured in the laboratory

% propagation of uncertainty from gain and offset to the calibrated
% radiance: for details see textbook on uncertainty
% the uncertainty is an absolute uncertainty
uL_calibrated = (uDN_.^2./gain.^2 + u_rad_coeffs.u_offset'.^2./gain.^2 + L_calibrated.^2.*u_rad_coeffs.u_gain'.^2./gain.^2 + 2*u_rad_coeffs.u_gain_offset.*L_calibrated./gain.^2).^(0.5);


fh=figure;
errorbar(wvl_instr, L_calibrated, uL_calibrated);
hold
plot(wvl_instr, L_calibrated, 'r', 'LineWidth', 3);
legend(escape_('uL_calibrated'),escape_('L_calibrated'))
xlabel('Wavelength [nm]');
ylabel('L [W/m^2/sr/nm]')
title_str = 'Calibrated Radiance with Uncertainty Envelope';
title(title_str)
print_pdf(fh, plot_dir, title_str);

% calculate relative uncertainty
uL_calibrated_rel = uL_calibrated ./ L_calibrated * 100;

fh=figure;
plot(wvl_instr, uL_calibrated_rel, 'r', 'LineWidth', 3);
xlabel('Wavelength [nm]');
ylabel('uL %]')
title_str = 'Relative Uncertainty of Calibrated Radiance at k=1';
title(title_str)
print_pdf(fh, plot_dir, title_str);







