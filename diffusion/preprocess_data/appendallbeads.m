%append all the c.mat files in a bigger structure
i=0; k=1;
while k==1
flag=menu('Select more beads','Yes','No');
if  flag==2

    if i==0
        k=0;
    else
        save('allbeads.mat', 'b');
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