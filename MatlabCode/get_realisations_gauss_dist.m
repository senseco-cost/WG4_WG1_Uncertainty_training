% creates realisations, either using a new random generated number array or
% using the supplied random number array
% Attention: sigma must be in absolute numbers and not a relative
% uncertainty
function [realisations, random_numbers] = get_realisations_gauss_dist(numSamples, mu, sigma, random_numbers)

    if nargin == 3
        random_numbers = randn(numSamples, length(mu));   
    end
    
         
        
    for i=1:numSamples
        
        realisations(:,i) = mu + sigma.*random_numbers(i,:)';
        
    end
    
    
    % figure
    % plot(mu + sigma.*random_numbers(i,:)')
    % hold
    % plot(mu, 'r');


end
