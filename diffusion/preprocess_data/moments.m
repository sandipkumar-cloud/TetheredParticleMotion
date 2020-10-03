% calculate the moments of distribution of all beads
% input is the overall rho of all beads
clear e
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);
delete([pathname,'momentallbeadsALL.xlsx']);
xlsfile = [pathname 'momentallbeadsALL.xlsx'];

m=0;
for j=1:length(b)
    if isempty(b(1,j).OK.indices)==0      
     for k=1:size(b(1,j).OK.rho_overall,2)
          m=m+1;
          for l=1:size(b(1,j).OK.rho_overall,1)
          e.ALL.rho_overall(l,m)=b(1,j).OK.rho_overall(l,k);
          end
         e.ALL.length(m)=size(b(1,j).OK.rho_overall,1);
         e.ALL.bead_number(m)=m;
     end
    end
end

for j=1:length(b)
    if isempty(b(1,j).REC.indices)==0
     for k=1:size(b(1,j).REC.rho_overall,2)
          m=m+1;
          for l=1:b(1,j).REC.length(k)
          e.ALL.rho_overall(l,m)=b(1,j).REC.rho_overall(l,k);
          end
         e.ALL.length(m)=b(1,j).REC.length(k);  
         e.ALL.bead_number(m)=m;
     end
    end
end


% rho and the length is saved for plotting rho or doing other comparisons
% of all the rho of one dna tether length in one buffer condition

e.xout=[-20:5:1000];
for i=1:length(e.ALL.length)
    y=e.ALL.rho_overall(1:e.ALL.length(i),i);
    n=hist(y,e.xout);
    for j=1:length(n)
        e.n_norm(j,i)=n(j)/sum(n);
    end
    e.mean(i)=mean(e.ALL.rho_overall(1:e.ALL.length(i),i));
    e.var(i)=var(e.ALL.rho_overall(1:e.ALL.length(i),i));
    e.skewness(i)=skewness(e.ALL.rho_overall(1:e.ALL.length(i),i));
    e.kurtosis(i)=kurtosis(e.ALL.rho_overall(1:e.ALL.length(i),i));
        e.mean_moment(1,i)=mean(e.ALL.rho_overall(1:e.ALL.length(i),i));
    e.mean_moment(2,i)=var(e.ALL.rho_overall(1:e.ALL.length(i),i));
    e.mean_moment(3,i)=skewness(e.ALL.rho_overall(1:e.ALL.length(i),i));
    e.mean_moment(4,i)=kurtosis(e.ALL.rho_overall(1:e.ALL.length(i),i));
    
end

xlswrite(xlsfile,e.mean,'mean','A2');
xlswrite(xlsfile,e.var,'variance','A2');
xlswrite(xlsfile,e.skewness,'skewness','A2');
xlswrite(xlsfile,e.kurtosis,'kurtosis','A2');
xlswrite(xlsfile,e.mean_moment,'central moments','A2');
xlswrite(xlsfile,e.xout,'xout','A1');
xlswrite(xlsfile,e.n_norm,'n_norm','A2');

xlswrite(xlsfile,e.ALL.bead_number,'n_norm','A1');
xlswrite(xlsfile,e.ALL.bead_number,'mean','A1');
xlswrite(xlsfile,e.ALL.bead_number,'variance','A1');
xlswrite(xlsfile,e.ALL.bead_number,'skewness','A1');
xlswrite(xlsfile,e.ALL.bead_number,'kurtosis','A1');
xlswrite(xlsfile,e.ALL.bead_number,'central moments','A1');

fig=figure(1);
plot(e.xout,e.n_norm);
saveas(fig,'fig_all_rho')
save('momentsALL.mat', 'e');

