% selected beads
clear all
close all
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);
delete([pathname,'rho_distribution_selected.xlsx']);
xlsfile = [pathname 'rho_distribution_selected.xlsx']; % creation of excel file for rho distribution

for m=1:length(d.selected_bead_number)
    n=d.selected_bead_number(m);
    d.n_rho_norm_selected(m,:)=d.n_rho_norm(n,:);
    d.selected_length(m)=d.length(n);
end

xlswrite(xlsfile,d.rho_bins,'rho_dist','B1');
xlswrite(xlsfile,d.n_rho_norm_selected,'rho_dist','B2');
xlswrite(xlsfile,d.selected_bead_number,'rho_dist','A2');

fig=figure(1);
plot(d.rho_bins,d.n_rho_norm_selected);
saveas(fig,'fig_selected_rho')
save('analyzed8s.mat', 'd');
clear all
close all