% calculate the moments of distribution of all beads
% input is the overall rho of all beads
clear e f h
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name2=[pathname, filename];
load(filename);
delete([pathname,'momentallbeadsGOOD.xlsx']);
xlsfile = [pathname 'momentallbeadsGOOD.xlsx'];

for m=1:length(f)
    n=f(m);
    
    h.GOOD.rho_overall(:,m)=e.ALL.rho_overall(:,n);
    h.GOOD.x(:,m)=e.ALL.x(:,n);    
    h.GOOD.y(:,m)=e.ALL.y(:,n);
%     h.GOOD.time(:,m)=e.ALL.time(:,n);   
    h.GOOD.length(m)=e.ALL.length(n);
     h.GOOD.status(m)=e.ALL.status(n);
    h.GOOD.bead_number(m)=n;
end


% rho and the length is saved for plotting rho or doing other comparisons
% of all the rho of one dna tether length in one buffer condition
clear n
h.xout=[-20:5:1000];
for i=1:length(h.GOOD.length)
    y=h.GOOD.rho_overall(1:h.GOOD.length(i),i);
    n=hist(y,h.xout);
    for j=1:length(n)
        h.n_norm(j,i)=n(j)/sum(n);
    end
    h.mean(i)=mean(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.var(i)=var(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.skewness(i)=skewness(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.kurtosis(i)=kurtosis(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.mean_moment(1,i)=mean(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.mean_moment(2,i)=var(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.mean_moment(3,i)=skewness(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    h.mean_moment(4,i)=kurtosis(h.GOOD.rho_overall(1:h.GOOD.length(i),i));
    
end

xlswrite(xlsfile,h.mean,'mean','A2');
xlswrite(xlsfile,h.var,'variance','A2');
xlswrite(xlsfile,h.skewness,'skewness','A2');
xlswrite(xlsfile,h.kurtosis,'kurtosis','A2');
xlswrite(xlsfile,h.mean_moment,'central moments','A2');
xlswrite(xlsfile,h.xout,'xout','A1');
xlswrite(xlsfile,h.n_norm,'n_norm','A2');

xlswrite(xlsfile,h.GOOD.bead_number,'n_norm','A1');
xlswrite(xlsfile,h.GOOD.bead_number,'mean','A1');
xlswrite(xlsfile,h.GOOD.bead_number,'variance','A1');
xlswrite(xlsfile,h.GOOD.bead_number,'skewness','A1');
xlswrite(xlsfile,h.GOOD.bead_number,'kurtosis','A1');
xlswrite(xlsfile,h.GOOD.bead_number,'central moments','A1');

fig=figure(1);
plot(h.xout,h.n_norm);
saveas(fig,'fig_good_rho')
save('momentsGOOD.mat', 'h');

