%create individual directories.
path=('E:\Data2\Ely_May28th\Output charts');
cd(path);
d=dir();

for i=73:81
    for j=1:7
        for k=1:70
            
    A=importdata(d(i).name);
end