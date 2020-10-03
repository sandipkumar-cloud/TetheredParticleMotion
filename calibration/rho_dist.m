% rho distribution of beads and save the distribution for clustering

% chose file to analyze (choose analyzed8s.mat)
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
addpath(pathname)
name=[pathname, filename];
load(filename);
delete([pathname,'rho_distribution.xlsx']);
xlsfile = [pathname 'rho_distribution.xlsx']; % creation of excel file for rho distribution

d.rho_bins=0:5:1000; % creation of bins

for i=1:length(d.length)
    y=d.rho(1:d.length(i),i);
    n=hist(y,d.rho_bins);
    for j=1:length(n)
        d.n_rho_norm(i,j)=n(j)/sum(n);
    end  
end

d.bead_number=1:length(d.length);

xlswrite(xlsfile,d.rho_bins,'rho_dist','B1');
xlswrite(xlsfile,d.n_rho_norm,'rho_dist','B2');
xlswrite(xlsfile,d.bead_number.','rho_dist','A2');

fig=figure(1);
plot(d.rho_bins,d.n_rho_norm);
saveas(fig,'fig_all_rho')
save('analyzed8s.mat', 'd');
clear all
close all
