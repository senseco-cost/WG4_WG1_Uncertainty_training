%% Monte Carlo Example



% draw from a Gaussian distribution

N = 100;
mu = 42;
sigma = 1;

rel_sigma = sigma / mu * 100;

[realisations, random_numbers] = get_realisations_gauss_dist(N, mu, sigma);


fh=figure;
plot(1:N, ones(N,1)*mu, 'r');
hold
plot(realisations, 'o')
title_str = ['Realisations of 42 with uncertainty 1, N = ' num2str(N)];
title(title_str);
xlabel('# Realisations')
ylabel('Realisations of the answer to the universe, life and everyting')

print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);


fh=figure;
histogram(realisations)
title_str = ['Histogram of Realisations of 42 with uncertainty 1, N = ' num2str(N)];
title(title_str)

print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);



N = 1000;

[realisations, random_numbers] = get_realisations_gauss_dist(N, mu, sigma);


fh=figure;
plot(1:N, ones(N,1)*mu, 'r');
hold
plot(realisations, 'o')
title_str = ['Realisations of 42 with uncertainty 1, N = ' num2str(N)];
title(title_str);
xlabel('# Realisations')
ylabel('Realisations of the answer to the universe, life and everyting')

print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);


fh=figure;
histogram(realisations)
title_str = ['Histogram of Realisations of 42 with uncertainty 1, N = ' num2str(N)];
title(title_str)

print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);


N = 10000;

[realisations, random_numbers] = get_realisations_gauss_dist(N, mu, sigma);


fh=figure;
plot(1:N, ones(N,1)*mu, 'r');
hold
plot(realisations, 'o')
title_str = ['Realisations of 42 with uncertainty 1, N = ' num2str(N)];
title(title_str);
xlabel('# Realisations')
ylabel('Realisations of the answer to the universe, life and everyting')

print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);


fh=figure;
histogram(realisations)
title_str = ['Histogram of Realisations of 42 with uncertainty 1, N = ' num2str(N)];
title(title_str)

print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);



% uncertainty propagation



% simple model to draw from the generated parameter values

output = zeros(size(realisations)); % allocate output vector

% loop over all realisations and calculate the model output
for i=1:length(realisations)
    
    parameter = realisations(i);
    
    output(i) = (log(parameter))^2;
    
end


linewidth = 2;
FontSize = 12;
TitleFontSize = 15;

% plot output as vector
fh=figure;
plot(output, 'o')
xlabel('Count', 'FontSize', FontSize);
ylabel('Model Value', 'FontSize', FontSize);
title_str = 'Values produced by Model for Parameter Realisations';
title(title_str, 'FontSize', TitleFontSize-1);
print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);

% plot output as histogram
fh=figure;
hist(output, 25)
xlabel('Value', 'FontSize', FontSize);
ylabel('Frequency', 'FontSize', FontSize);
title_str = 'Histogram of Model Output';
title(title_str, 'FontSize', TitleFontSize);
print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', title_str);



% N = 100000;
% 
% [realisations, random_numbers] = get_realisations_gauss_dist(N, mu, sigma);
% output = zeros(size(realisations)); % allocate output vector
% 
% % loop over all realisations and calculate the model output
% for i=1:length(realisations)
%     
%     parameter = realisations(i);
%     
%     output(i) = (log(parameter))^2;
%     
% end
% 

% get statistics of the output
mean_output = mean(output)
stddev_output = std(output)

rel_sigma_output = stddev_output / mean_output * 100;




%% analytical deduction

x=min(realisations):max(realisations);
y = (log(x)).^2;


fh=createfigureOfLinearFit(x, y);
print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', 'Linearity of log squared');


x=20:50;
y = (log(x)).^2;


fh=createfigureOfLinearFit(x, y);
print_pdf(fh, '/Users/andyhueni/Data/Studies/RSL/COST/SENSECO/training school 2019/Plots/', 'Linearity of log squared - non linearity example');


c = (2*log(mu)/mu);

u2_y = c^2 * sigma^2;
u_y = u2_y ^ (0.5)

delta = stddev_output - u_y




