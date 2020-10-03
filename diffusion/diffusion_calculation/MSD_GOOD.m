%analyse all beads
clear all
[filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
name=[pathname, filename];
load(filename);
v=h.GOOD;
wmax=10;
a=1:wmax;
b=2*2.^(a-1);
v.t=20.*b;%time window in milliseconds

for j=1:length(v.length)
    x_even=v.x(2:2:(2*abs(v.length(j)/2)),j);
    y_even=v.y(2:2:(2*abs(v.length(j)/2)),j);
    v.mean_x_even(j)=mean(x_even(1:abs(v.length(j)/2)));
    v.mean_y_even(j)=mean(y_even(1:abs(v.length(j)/2)));
    x_odd=v.x(1:2:(2*abs(v.length(j)/2)-1),j);
    y_odd=v.y(1:2:(2*abs(v.length(j)/2)-1),j);
    v.mean_x_odd(j)=mean(x_odd(1:abs(v.length(j)/2)));
    v.mean_y_odd(j)=mean(y_odd(1:abs(v.length(j)/2)));
    clear x_even y_even x_odd y_odd
    for n=1:wmax
        f=2*2^(n-1);
        %             even frames
        m=1;
        for l=2*f:2:(2*abs(v.length(j)/2))
            v.sqDx_even(j,n,m)=(v.x(l,j)-v.x((l-f),j))^2;
            v.sqDy_even(j,n,m)=(v.y(l,j)-v.y((l-f),j))^2;
            v.sqD_even(j,n,m)=v.sqDx_even(j,n,m)+v.sqDy_even(j,n,m);
            v.x_even(j,n,m)=v.x(l,j)-v.mean_x_even(j);
            v.y_even(j,n,m)=v.y(l,j)-v.mean_y_even(j);
            v.rho_even(j,n,m)=(v.x_even(j,n,m)^2+v.y_even(j,n,m)^2)^(1/2);
            m=m+1;
        end
        %             odd frames
        o=1;
        for l=2*f-1:2:(2*abs(v.length(j)/2)-1)
            v.sqDx_odd(j,n,o)=(v.x(l,j)-v.x((l-f),j))^2;
            v.sqDy_odd(j,n,o)=(v.y(l,j)-v.y((l-f),j))^2;
            v.sqD_odd(j,n,o)=v.sqDx_odd(j,n,o)+v.sqDy_odd(j,n,o);
            v.x_odd(j,n,o)=v.x(l,j)-v.mean_x_odd(j);
            v.y_odd(j,n,o)=v.y(l,j)-v.mean_y_odd(j);
            v.rho_odd(j,n,o)=(v.x_odd(j,n,o)^2+v.y_odd(j,n,o)^2)^(1/2);
            o=o+1;
        end
    end
end


% calculation of MSD of each bead at each time intervals
for j=1:length(v.length)
    
    for n=1:wmax
        f=2*2^(n-1);
        %             even frames
        
        v.MSD_x_even(j,n)=mean(v.sqDx_even(j,n,1:(abs(v.length(j)/2)-f+1)),3);
        v.MSD_y_even(j,n)=mean(v.sqDy_even(j,n,1:(abs(v.length(j)/2)-f+1)),3);
        v.MSD_even(j,n)=v.MSD_x_even(j,n)+v.MSD_y_even(j,n);
        %             odd frames
        
        v.MSD_x_odd(j,n)=mean(v.sqDx_odd(j,n,1:(abs(v.length(j)/2)-f+1)),3);
        v.MSD_y_odd(j,n)=mean(v.sqDy_odd(j,n,1:(abs(v.length(j)/2)-f+1)),3);
        v.MSD_odd(j,n)=v.MSD_x_odd(j,n)+v.MSD_y_odd(j,n);
        
    end
end


save('MSD_GOOD.mat','v');
clear all

