%% RTM_and_uncertainty propagation_Ex_1_solved.m
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
% % Download from http://teledetection.ipgp.jussieu.fr/prosail/
addpath('../PROSAIL_D_MATLAB_2017')
%% Load data
[num,txt] = xlsread(['Ex1_TableLeafParam.csv']);
% % Assign data to variables
for i=1:length(txt)
    eval(sprintf('%s = [%f,%f,%f];',txt{i},num(:,i)));
end
clear txt num

%% Define functions to compute leaf radious, perimeter and area
% % Assume that the largest radious is twice the minor one in the eliptical
% leaves
% % Minor radious as a function of the leaf Area
rf = @(A)sqrt(A./(2*pi)); 
% % Elipse area as a function of the minor radious of the elipse
Af = @(x)pi.*x.*(2.*x);
% % Elipse perimeter as a function of the minor radious of the elipse
Df = @(x)2*pi*sqrt((x.^2 + (2.*x).^2)/2);
% % Function that converts centimeters to pixels
cm2pix = @(x,res)(res.*x)./2.54;

% % Estimate the minor radious of each eliptical leaf
r = rf(A);

%% Define uncertainties (see .ppt) and scanner resolution
res = 300; % dpi; divide DPI by 2.54 cm/inch to get pixel per cm
ucab = 3;
ucar = 1.5;
uw = 0.02;
% % Uncertainty of area due to indetermination of border pixels
uA = A.*( Df(cm2pix(r,res))./Af(cm2pix(r,res)) );

%% Uncertainty of water and dry matter content
Cm = wd./A;
ucm = sqrt((uw./wd).^2 + (uA./A).^2);

deltaw = (wf - wd);
udeltaw = sqrt(uw.^2 + uw.^2);
Cw = deltaw./A;
ucw = sqrt((udeltaw./deltaw).^2 + (uA./A).^2);

%% Simulate leaf optical properties
% % Number of simulations per leaf
n_mc = 100;
% % Define wavelength
wl = [400:2500]';
for i=1:length(N) % Number of leaves
	% % Preallocates 
	[rho,tau] = deal(nan(length(wl),n_mc));
    
	% %  Draw distributions of the parameters
    Cabi = normrnd_truncated(Cab(i),ucab,[n_mc],0,100);
    Cari = normrnd_truncated(Car(i),ucar,[n_mc],0,40);
       
    if true        
        Cmi = normrnd_truncated(Cm(i),Cm(i).*ucm(i),[n_mc],0,0.03);
        Cwi = normrnd_truncated(Cw(i),Cw(i).*ucw(i),[n_mc],0,0.04);
    else
        % % Alternativey Monte Carlo can be applied to the weights and the
        % % areas
        wfi = normrnd_truncated(wf(i),uw,[n_mc],0,inf);
        wdi = normrnd_truncated(wd(i),uw,[n_mc],0,inf);

        Ai = normrnd_truncated(A(i),uA(i),[n_mc],0,inf);

        % % Compute remaining leaf parameters
        Cmi = wdi./Ai;
        Cwi = (wfi - wdi)./Ai;   
    end
    
% % Run the model for each simulation
   for j = 1:n_mc
       LRT = prospect_DB(N(i),Cabi(j),Cari(j),0,0,Cwi(j),Cmi(j));
       rho(:,j) = LRT(:,2);
       tau(:,j) = LRT(:,3);
   end
   
% % Plot the result
   clf; set(gcf,'Color','w'); hold on; grid on; axis([400 2500 0 1])
   p1 = plot(wl,rho,'Color',[.7 .7 .7]);
   plot(wl,1-tau,'Color',[.7 .7 .7]);
   p2 = plot(wl,nanmean(rho,2),'Color',[.2 .8 .2]);
   p3 = plot(wl,nanmean(rho,2) - nanstd(rho,[],2),'--','Color',[.2 .2 .8]);   
   plot(wl,nanmean(rho,2) + nanstd(rho,[],2),'--','Color',[.2 .2 .8]);
   plot(wl,1-nanmean(tau,2),'Color',[.2 .8 .2]);
   plot(wl,1-nanmean(tau,2) - nanstd(tau,[],2),'--','Color',[.2 .2 .8]);   
   plot(wl,1-nanmean(tau,2) + nanstd(tau,[],2),'--','Color',[.2 .2 .8]);
   title(sprintf('Leaf %d',i));
   xlabel('Wavelength (nm)'); ylabel('\it\rho\rm or 1-\it\tau\rm')
   
   l = legend([p1(1),p2(1),p3(1)],'Sim.','Mean','\itu\rm',...
       'location','best');
   saveas(figure(1),[pwd,'/Done_PROSPECT_sim_',num2str(i),'.png'])
   pause(1)
end

%% Remove paths
rmpath('../PROSAIL_D_MATLAB_2017')
