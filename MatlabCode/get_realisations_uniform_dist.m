% creates realisations, either using a new random generated number array or
% using the supplied random number array
% Samples are taken from a uniform distribution
function [realisations, random_numbers] = get_realisations_uniform_dist(numSamples, mu, interval, random_numbers)

    if nargin == 3
        random_numbers = rand(numSamples, length(mu));   
    end
    
         
        
    for i=1:numSamples
        
        realisations(:,i) =  (interval.*(random_numbers(i,:)) +(mu-interval/2))'; % centre on zero, then add mu to centre on actual mean
        
    end
    
    do_plot = true;
    if do_plot
        
%         figure
%         plot(squeeze(interval.*(random_numbers)))
%         
        figure
        plot(realisations)
        hold
        plot(mu, 'r');
        
        figure
        histogram(realisations);
        
        min(realisations);
        max(realisations)
    end

end
