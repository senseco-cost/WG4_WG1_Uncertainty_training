%% RTM_and_uncertainty propagation_Ex_2_solved.m
% %     This is script contains a solution proposed for the Exercise 1 of 
% % the Session 2: "Uncertainty and variability propagation trough RTMs" of 
% % the WG4 and WG1 training school "Traceability Chains, Uncertainty 
% % Propagation and Calibration/Validation" organized by COST Action 
% % CA17134 - Optical synergies for spatiotemporal SENSing of scalable 
% % ECOphysiological traits (SENSECO) between the 18-22 November 2019 in
% % CETAL, Magurele, Romania.
% % 
% % Author: 
% %     Javier Pacheco-Labrador, PhD  (jpacheco@bgc-jena.mpg.de)
% %     MPI-BGC, Jena, Germany
% % 

close all
clear all
clc
addpath('../PROSAIL_D_MATLAB_2017')

%% Load data
% % Pedefine the number of samples
n_subplots = 16;
% % Load leaf metrics
[num,txt] = xlsread(['Ex2_TableLeafParam.csv']);
% % Assign data to variables
for i=1:length(txt)
    eval(sprintf(['%s = [%f',repmat(',%f',1,n_subplots-1),'];'],txt{i},num(:,i)));
end
clear txt num


% % Load soil refletance factor
rsoil = csvread('Soil.csv'); % clf; plot(rsoil(:,1),rsoil(:,2));
% % Remove wavelength
rsoil = rsoil(:,2);


% % Load data for Leaf Area Index
[num,txt] = xlsread(['Ex2_TableLAI.csv']);
% % Assign data to variables
for i=1:length(txt)
    eval(sprintf(['%s = [%f',repmat(',%f',1,n_subplots-1),'];'],txt{i},num(:,i)));
end
clear txt num


% % Load Bottom of the Atmosphere irradiances
[num,txt] = xlsread(['Ex2_BOAirradiance.csv']);
% % Assign data to variables
for i=1:length(txt)
    eval(sprintf(['%s = transpose([%f',repmat(',%f',1,length(rsoil)-1),']);'],txt{i},num(:,i)));
end
clear txt num


%% Define canopy and sun-view parameters 
tts = 30;
tto = 0;
psi = 0;

LIDFa = -0.35;
LIDFb = -0.15;

hot_spot = 0.01;
Cbrown = 0;
Cant= 0;


%% Define ancillary functions
% % Assume that the largest radious is twice the minor one in the eliptical
% leaves
% % Minor radious as a function of the leaf Area
rf = @(A)sqrt(A./(2*pi)); 
% % Elipse perimeter as a function of the minor radious of the elipse
Af = @(x)pi.*x.*(2.*x);
% % Elipse area  as a function of the minor radious of the elipse
Df = @(x)2*pi*sqrt((x.^2 + (2.*x).^2)/2);
% % Function that converts centimeters to pixels
cm2pix = @(x,res)(res.*x)./2.54;

% % Estimate the minor radious of each eliptical leaf
r = rf(A);

%% Define LAI model and functions
% % Define aPAR
aPARf = @(Ein,Etran)(Ein./Etran);
% % Define model predicting LAI from aPAR
LAIf = @(aPAR)(-1.4829.*log(aPAR) + 0.0246);

%% Define uncertainties (see .ppt) and scanner resolution
res = 300; % dpi; divide DPI by 2.54 cm/inch to get pixel per cm
ucab = 3;
ucar = 1.5;
uw = 0.01;
% % Uncertainty of area due to indetermination of border pixels
uA = A.*( Df(cm2pix(r,res))./Af(cm2pix(r,res)) );

% % Define the spatial variability of the parameters in relative terms
sigma_spatial = 0.05;

% % Define relative uncertainty of the irradiance sensor for aPAR
uE = 0.05;
%% Simulate leaf optical properties
% % Number of simulations per leaf
n_mc = 100;
% % Define wavelength
wl = [400:2500]';
for i=1:n_subplots % Number of leaves
    fprintf('Subplot %d',i)
	% % Preallocates 
	[rho,tau,Rtoc] = deal(nan(length(wl),n_mc));
    
% % Draw distributions of the leaf parameters.
% % Here the measurement uncertainty with the spatial variability of the
% % variable are combined assuming no correlation (sum of squares)
% %     Consider that: norm = sqrt(sum(x.^2))
    Cabi = normrnd_truncated(Cab(i),norm([ucab,Cab(i).*sigma_spatial]),[n_mc],0,100);
    Cari = normrnd_truncated(Car(i),norm([ucar,Car(i).*sigma_spatial]),[n_mc],0,40);
   
    wfi = normrnd_truncated(wf(i),norm([uw,wf(i).*sigma_spatial]),[n_mc],0,inf);
    wdi = normrnd_truncated(wd(i),norm([uw,wd(i).*sigma_spatial]),[n_mc],0,inf);
   
    Ai = normrnd_truncated(A(i),norm([uA(i),A(i).*sigma_spatial]),[n_mc],0,inf);
   
% % Compute remaining leaf parameters
    Cmi = wdi./Ai;
    Cwi = (wfi - wdi)./Ai; 
    
% % Estimate LAI 
% % First draw the aPAR distribution
    aPAR = aPARf(Ein(i).*normrnd(1,.05,[n_mc,1]),Etran(i).*normrnd(1,.05,[n_mc,1]));
% % Then predict LAI and add uncertainty of the preductive model
    LAIi = LAIf(aPAR) + normrnd(1,.2522,[n_mc,1]);
    
% % Run the model for each simulation
   for j = 1:n_mc
% %    Run PROSPECT separately for the sake of getting the plots
       LRT = prospect_DB(N(i),Cabi(j),Cari(j),0,0,Cwi(j),Cmi(j));
       rho(:,j) = LRT(:,2);
       tau(:,j) = LRT(:,3);
       
% %    Run PROSAIL
       [rdot,rsot]=PRO4SAIL(N(i),Cabi(j),Cari(j),Cant,Cbrown,Cwi(j),Cmi(j),...
           LIDFa,LIDFb,1,LAIi(j),hot_spot,tts,tto,psi,rsoil);    
       
% %     Compute the observed reflectance factor according to the fractions
% %     of direct and diffuse light
        Rtoc(:,j) = (rdot.*Edif + rsot.*Edir)./(Edir + Edif);
   end
   
% % Plot the result
   figure(1); clf; set(gcf,'Color','w'); hold on; grid on; axis([400 2500 0 1])
   p1 = plot(wl,rho,'Color',[.7 .7 .7]);
   plot(wl,1-tau,'Color',[.7 .7 .7]);
   p2 = plot(wl,nanmean(rho,2),'Color',[.2 .8 .2]);
   p3 = plot(wl,nanmean(rho,2) - nanstd(rho,[],2),'--','Color',[.2 .2 .8]);   
   plot(wl,nanmean(rho,2) + nanstd(rho,[],2),'--','Color',[.2 .2 .8]);
   plot(wl,1-nanmean(tau,2),'Color',[.2 .8 .2]);
   plot(wl,1-nanmean(tau,2) - nanstd(tau,[],2),'--','Color',[.2 .2 .8]);   
   plot(wl,1-nanmean(tau,2) + nanstd(tau,[],2),'--','Color',[.2 .2 .8]);
   title(sprintf('Leaf subplot %d',i));
   xlabel('Wavelength (nm)'); ylabel('\it\rho\rm or 1-\it\tau\rm')
   
   l = legend([p1(1),p2(1),p3(1)],'Sim.','Mean','\itu\rm',...
       'location','best');
   saveas(figure(1),[pwd,'/Ex2_PROSPECT_sim_',num2str(i),'.jpg'])
   
   
   figure(2);  clf; set(gcf,'Color','w'); hold on; grid on; axis([400 2500 0 .5])
   p1 = plot(wl,Rtoc,'Color',[.7 .7 .7]);
   p2 = plot(wl,nanmean(Rtoc,2),'Color',[.2 .8 .2]);
   p3 = plot(wl,nanmean(Rtoc,2) - nanstd(Rtoc,[],2),'--','Color',[.2 .2 .8]);
   plot(wl,nanmean(Rtoc,2) + nanstd(Rtoc,[],2),'--','Color',[.2 .2 .8]);
   title(sprintf('Canopy subplot %d',i));
   xlabel('Wavelength (nm)'); ylabel('\itR\rm')
   
   l = legend([p1(1),p2(1),p3(1)],'Sim.','Mean','\itu\rm',...
       'location','best');
   saveas(figure(2),[pwd,'/Ex2_PROSAIL_sim_',num2str(i),'.jpg'])
   
end


%% Remove paths
rmpath('../PROSAIL_D_MATLAB_2017')

