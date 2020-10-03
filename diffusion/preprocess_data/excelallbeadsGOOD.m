%analyse all beads
clear d
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name1=[pathname, filename];
load(filename);
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name2=[pathname, filename];
load(filename);
delete([pathname,'excelallbeadsGOOD.xlsx']);
xlsfile = [pathname 'excelallbeadsGOOD.xlsx'];
for m=1:length(f)
    n=f(m);
    for l=1:length(d.ALL.filter_rect)
    g.GOOD.rho_avg_avg(m,l)=d.ALL.rho_avg_avg(n,l);
    g.GOOD.rho_rms_avg(m,l)=d.ALL.rho_rms_avg(n,l);
    g.GOOD.rho_avg_std(m,l)=d.ALL.rho_avg_std(n,l);
    g.GOOD.rho_rms_std(m,l)=d.ALL.rho_rms_std(n,l);
    g.GOOD.rho_stdt(m,l)=d.ALL.rho_stdt(n,l);
    g.GOOD.rho_median_median(m,l)=d.ALL.rho_median_median(n,l);
    g.GOOD.rho_median_mad(m,l)=d.ALL.rho_median_mad(n,l);
    end
    g.GOOD.bead_number(m,1)=m;
    
end


g.GOOD.filter_rect=d.ALL.filter_rect;
g.GOOD.rho_median_all=median(g.GOOD.rho_median_median);
g.GOOD.rho_mad_all=mad(g.GOOD.rho_median_median);

xlswrite(xlsfile,g.GOOD.rho_avg_avg,'rho avg avg','B3');
xlswrite(xlsfile,g.GOOD.rho_rms_avg,'rho rms avg','B3');
xlswrite(xlsfile,g.GOOD.rho_avg_std,'rho avg std','B3');
xlswrite(xlsfile,g.GOOD.rho_rms_std,'rho rms std','B3');
xlswrite(xlsfile,g.GOOD.rho_stdt,'rho stdt','B3');
xlswrite(xlsfile,g.GOOD.rho_median_median,'rho median median','B3');
xlswrite(xlsfile,g.GOOD.rho_median_mad,'rho median mad','B3');
xlswrite(xlsfile,g.GOOD.rho_median_all,'rest','A2');
xlswrite(xlsfile,g.GOOD.rho_mad_all,'rest','A3');
xlswrite(xlsfile,g.GOOD.filter_rect,'rest','A1');

xlswrite(xlsfile,g.GOOD.bead_number,'rho avg avg','A3');
xlswrite(xlsfile,g.GOOD.bead_number,'rho rms avg','A3');
xlswrite(xlsfile,g.GOOD.bead_number,'rho avg std','A3');
xlswrite(xlsfile,g.GOOD.bead_number,'rho rms std','A3');
xlswrite(xlsfile,g.GOOD.bead_number,'rho stdt','A3');
xlswrite(xlsfile,g.GOOD.bead_number,'rho median median','A3');
xlswrite(xlsfile,g.GOOD.bead_number,'rho median mad','A3');

save('allbeadsanalyzedGOOD.mat', 'g');
clear all

    