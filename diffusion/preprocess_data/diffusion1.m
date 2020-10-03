%analyse all beads
clear v
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);
delete([pathname,'displacementALL.xlsx']);
xlsfile = [pathname 'displacementALL.xlsx'];
m=0;
for j=1:length(b)
    if isempty(b(1,j).OK.indices)==0
        for k=1:size(b(1,j).OK.rho_avg_avg,2)
            m=m+1;
            for n=1:size(b(1,j).OK.beadsx,1)
                v.ALL.beadsx(n,m)=b(1,j).OK.beadsx(n,k);
                v.ALL.beadsy(n,m)=b(1,j).OK.beadsy(n,k);
            end
            for n=1:(size(b(1,j).OK.beadsx,1)-1)
                v.ALL.displacement_square(n,m)=((v.ALL.beadsx(n,m)-v.ALL.beadsx((n+1),m))^2+(v.ALL.beadsy(n,m)-v.ALL.beadsy((n+1),m))^2);
                v.ALL.displacement(n,m)=sqrt((v.ALL.beadsx(n,m)-v.ALL.beadsx((n+1),m))^2+(v.ALL.beadsy(n,m)-v.ALL.beadsy((n+1),m))^2);
            end
            v.ALL.displacement_square_mean(m,1)=mean(v.ALL.displacement_square(1:(size(b(1,j).OK.beadsx,1)-1),m));
            v.ALL.diffusion(m,1)=v.ALL.displacement_square_mean(m,1)/(4*0.02);
            v.ALL.bead_number(m,1)=m;
        end
    end
end
v.ALL.diffusion_no_outliers=removeoutliers(v.ALL.diffusion);
xlswrite(xlsfile,v.ALL.diffusion_no_outliers,'diffusion','C3');



xlswrite(xlsfile,v.ALL.displacement_square_mean,'avg displacement square','B3');
xlswrite(xlsfile,v.ALL.bead_number,'avg displacement square','A3');
xlswrite(xlsfile,v.ALL.diffusion,'avg displacement square','C3');


save('displacementALL.mat','v');
clear all

