function X = normrnd_truncated(mu,sigma,tam,l,u)
% % HELP: This function samples form a normal truncated distribution. This 
% % function is desinged for vectors only.
% %     This is function was done for the   "Uncertainty and variability 
% % propagation trough RTMs" of the WG4 and WG1 training school "Traceability 
% % Chains, Uncertainty Propagation and Calibration/Validation" organized 
% % by COST Actio CA17134 - Optical synergies for spatiotemporal SENSing of 
% % scalable ECOphysiological traits (SENSECO) between the 18-22 November 
% % 2019 in CETAL, Magurele, Romania.
% % 
% % Author: 
% %     Javier Pacheco-Labrador, PhD  (jpacheco@bgc-jena.mpg.de)
% %     MPI-BGC, Jena, Germany


% % Preallocates
X = nan([tam,1]);

% % Generates distribution
i = 1; 
j = 1;
while i < tam
% % Draw a sample the double of the desired size
    xi = normrnd(mu,sigma,[tam*2,1]);
    
% % Keep values within the boundaries
    xi(xi<l | xi>u) = [];
% % Count the values within the boundaries
    lx = length(xi);
    
% % Find the last value that will use to fill up the table 
    j = min(tam,i + lx - 1);
% % Fill up the table with the valid values
    j = min(tam,i + lx - 1); 
    X(i:j) = xi(1:j-i+1);
% % Update the stating point of he next round
    i = j+1;
end
end

