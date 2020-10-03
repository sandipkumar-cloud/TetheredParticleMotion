%append all the c.mat files in a bigger structure
i=0; k=1;
while k==1
    flag=menu('Select more beads','Yes','No');
    if  flag==2
        
        if i==0
            k=0;
        else
            m=0;
            for j=1:length(b)
                if isempty(b(1,j).OK.indices)==0
                    for n=1:size(b(1,j).OK.beadsx,2)
                        m=m+1;
                        e.x(1:size(b(1,j).OK.beadsx,1),m)=b(1,j).OK.beadsx(1:size(b(1,j).OK.beadsx,1),n);
                        e.y(1:size(b(1,j).OK.beadsx,1),m)=b(1,j).OK.beadsy(1:size(b(1,j).OK.beadsx,1),n);
                        e.t(1:size(b(1,j).OK.beadsx,1),m)=b(1,j).OK.time(1,1:size(b(1,j).OK.beadsx,1));
                        e.length(m)=size(b(1,j).OK.beadsx,1);
                        e.SOURCE(m).path=strcat(b(1,j).SOURCE);
                    end
                end
            end
            
            for j=1:length(b)
                if isempty(b(1,j).REC.indices)==0
                    for n=1:size(b(1,j).REC.beadsx,2)
                        m=m+1;
                        e.x(1:b(1,j).REC.length(n),m)=b(1,j).REC.beadsx(1:b(1,j).REC.length(n),n);
                        e.y(1:b(1,j).REC.length(n),m)=b(1,j).REC.beadsy(1:b(1,j).REC.length(n),n);
                        e.t(1:b(1,j).REC.length(n),m)=b(1,j).REC.time(1:b(1,j).REC.length(n),n);
                        e.length(m)=b(1,j).REC.length(n);
                        e.SOURCE(m).path=strcat(b(1,j).SOURCE);
                    end
                end
            end
            
            
            save('xytlengthpath.mat', 'e');
            
            k=0;
        end
    end
    
    if flag==1
        i=i+1;
        [filename, pathname] = uigetfile({'*.mat';'*.*'},'File Selector');
        addpath(pathname)
        name=[pathname, filename];
        load(filename);
        b(i)=c;
        clear c filename pathname name
    end
end
clear all