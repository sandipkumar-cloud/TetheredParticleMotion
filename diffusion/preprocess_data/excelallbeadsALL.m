%analyse all beads
clear d
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);
delete([pathname,'excelallbeadsALL.xlsx']);
xlsfile = [pathname 'excelallbeadsALL.xlsx'];
m=0;
for j=1:length(b)
    if isempty(b(1,j).OK.indices)==0
     for k=1:size(b(1,j).OK.rho_avg_avg,2)
          m=m+1;
        for l=1:size(b(1,j).OK.rho_avg_avg,3)
            d.ALL.rho_avg_avg(m,l)=b(1,j).OK.rho_avg_avg(1,k,l);
            d.ALL.rho_rms_avg(m,l)=b(1,j).OK.rho_rms_avg(1,k,l);
            d.ALL.rho_avg_std(m,l)=b(1,j).OK.rho_avg_std(1,k,l);
            d.ALL.rho_rms_std(m,l)=b(1,j).OK.rho_rms_std(1,k,l);
            d.ALL.rho_stdt(m,l)=b(1,j).OK.rho_stdt(1,k,l);
            d.ALL.rho_median_median(m,l)=b(1,j).OK.rho_median_median(k,l);
            d.ALL.rho_median_mad(m,l)=b(1,j).OK.rho_median_mad(k,l);
            d.ALL.bead_number(m,1)=m;
        end   
     end
    end
end
for j=1:length(b)
    if isempty(b(1,j).REC.indices)==0
     for k=1:size(b(1,j).REC.rho_avg_avg,2)
          m=m+1;
        for l=1:size(b(1,j).REC.rho_avg_avg,3)
            d.ALL.rho_avg_avg(m,l)=b(1,j).REC.rho_avg_avg(1,k,l);
            d.ALL.rho_rms_avg(m,l)=b(1,j).REC.rho_rms_avg(1,k,l);
            d.ALL.rho_avg_std(m,l)=b(1,j).REC.rho_avg_std(1,k,l);
            d.ALL.rho_rms_std(m,l)=b(1,j).REC.rho_rms_std(1,k,l);
            d.ALL.rho_stdt(m,l)=b(1,j).REC.rho_stdt(1,k,l);
            d.ALL.rho_median_median(m,l)=b(1,j).REC.rho_median_median(k,l);
            d.ALL.rho_median_mad(m,l)=b(1,j).REC.rho_median_mad(k,l);
            d.ALL.bead_number(m,1)=m;
        end
             
     end
    end
end
filter_rect=[0.04,0.08,0.16,0.32,0.64,1.28,2.56,5.12,10.24,20.48,40.96,81.92];
d.ALL.filter_rect=filter_rect;
d.ALL.rho_median_all=median(d.ALL.rho_median_median);
d.ALL.rho_mad_all=mad(d.ALL.rho_median_median);

xlswrite(xlsfile,d.ALL.rho_avg_avg,'rho avg avg','B3');
xlswrite(xlsfile,d.ALL.rho_rms_avg,'rho rms avg','B3');
xlswrite(xlsfile,d.ALL.rho_avg_std,'rho avg std','B3');
xlswrite(xlsfile,d.ALL.rho_rms_std,'rho rms std','B3');
xlswrite(xlsfile,d.ALL.rho_stdt,'rho stdt','B3');
xlswrite(xlsfile,d.ALL.rho_median_median,'rho median median','B3');
xlswrite(xlsfile,d.ALL.rho_median_mad,'rho median mad','B3');
xlswrite(xlsfile,d.ALL.rho_median_all,'rest','A2');
xlswrite(xlsfile,d.ALL.rho_mad_all,'rest','A3');
xlswrite(xlsfile,d.ALL.filter_rect,'rest','A1');

xlswrite(xlsfile,d.ALL.bead_number,'rho avg avg','A3');
xlswrite(xlsfile,d.ALL.bead_number,'rho rms avg','A3');
xlswrite(xlsfile,d.ALL.bead_number,'rho avg std','A3');
xlswrite(xlsfile,d.ALL.bead_number,'rho rms std','A3');
xlswrite(xlsfile,d.ALL.bead_number,'rho stdt','A3');
xlswrite(xlsfile,d.ALL.bead_number,'rho median median','A3');
xlswrite(xlsfile,d.ALL.bead_number,'rho median mad','A3');

save('allbeadsanalyzedALL.mat', 'd');
clear all

    