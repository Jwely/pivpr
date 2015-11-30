%this script is a new approach at data viewing, 
%Opens all 10 relevant figures for a given station and combines then into a
%single condensed one page image without legends, just for profile
%comparison.

path='E:\Data2\Ely_May28th\Output charts\Fig Files';
cd(path);

%load up velocity profiles
d=dir('*.fig');
