
%this rather unelegant script is to test out some uncertainty propogation.
clc; clear all; format shortg

%==========================================================================
%                           1: User input values
%========================================================================== 
in = [-126.616, 96.8363, 0,...          %sample 12 sensitivity coefficients
    4.58828, 0.0591239, 3.08236,...     %These are of units (pixels/mm)
    -0.232408, 5.48947, 0.227232,...
    3.69093, -0.0735101, -3.30088,...
    0.252406, 4.97684, 0.247578];

% Constants                             indexes         eq number
% "X mm", "Y mm", "Z mm",               1  2  3         N/A
% "L dX/dx", "L dX/dy", "L dX/dz",      4  5  6         1
% "L dY/dx", "L dY/dy", "L dY/dz",      7  8  9         3
% "R dX/dx", "R dX/dy", "R dX/dz",      10 11 12        2
% "R dY/dx", "R dY/dy", "R dY/dz"       13 14 15        4

dpix=[4.5;-4.5;0.5;0.5];                % sample set of pixel displacements

fact = 50e-3;                           % time-space factor to convert all 
                                        % real displacements into units of 
                                        % (m/s)
                                        
var=[0.125^2, 0,0;...                   % diagonal uncertainty matrix, this  
    0, 0.125^2, 0;...                   % comes from 2nd order up-sampling
    0, 0, 0.125^2];                    

%==========================================================================
%         2:  Assembling coefficient and variance matrices
%========================================================================== 

%Matrix of coefficients             %equation numbers included in set
A=[in(4), in(5), in(6);...          %1 2 3
    in(10), in(11), in(12);...
    in(7), in(8), in(9)];

B=[in(4), in(5), in(6);...          %1 2 4
    in(10), in(11), in(12);...
    in(13), in(14), in(15)];

C=[in(4), in(5), in(6);...          %1 3 4
    in(7), in(8), in(9);...
    in(13), in(14), in(15)];

D=[in(10), in(11), in(12);...       %2 3 4
    in(7), in(8), in(9);...
    in(13), in(14), in(15)];

disp('------Variance-covariance matrices (m/s)------');
VarA = (A\var)*inv(A)'./fact            %expect high Y uncert
VarB = (B\var)*inv(B)'./fact            %expect high Y uncert
VarC = (C\var)*inv(C)'./fact            %expect high X uncert
VarD = (D\var)*inv(D)'./fact            %expect high X uncert

% The diagonal elements of each U matrix are the variances(uncertainties)
% for each set of equations. Note that matrices A and B corresponding to
% equation sets containing BOTH x pixel displacement equations, and have
% much lower uncertainties overall than C and D.

%==========================================================================
%    3:  Solving each set of equations for specified pixel displacements
%========================================================================== 

Sol_A = A\[dpix(1); dpix(2); dpix(3)]./fact ;
Sol_B = B\[dpix(1); dpix(2); dpix(4)]./fact ;
Sol_C = C\[dpix(1); dpix(3); dpix(4)]./fact ;
Sol_D = D\[dpix(2); dpix(3); dpix(4)]./fact ;

%==========================================================================
%    3:  Printing solutions and associated uncertainty values
%========================================================================== 

disp('------Solutions +/- stdev(sqrt variance)------'); format longg
fprintf('\n');

%uncertainty is in terms of standard deviation, so the square roots of the
%diagonal elements of the variance matrices are taken.

disp('Solution A (m/s)')
for i=1:3
fprintf('\t%5.6f \t+/- %5.6f \n',Sol_A(i),VarA(i,i)^0.5);
end
fprintf('\n');

disp('Solution B (m/s)')
for i=1:3
fprintf('\t%5.6f \t+/- %5.6f \n',Sol_B(i),VarB(i,i)^0.5);
end
fprintf('\n');

disp('Solution C (m/s)')
for i=1:3
fprintf('\t%5.6f \t+/- %5.6f \n',Sol_C(i),VarC(i,i)^0.5);
end
fprintf('\n');

disp('Solution D (m/s)')
for i=1:3
fprintf('\t%5.6f \t+/- %5.6f \n',Sol_D(i),VarD(i,i)^0.5);
end
















