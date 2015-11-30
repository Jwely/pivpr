%script to model potential flow within the ODU windtunnel
clc; clear all; close all;

H=0.914;        % tunnel height dimension
W=1.219;        % tunnel width dimension

% x portion of the profile

x=[0:0.001:W/2]';

for i=1:size(x,1)
        
        term(i,1)=1;
        term(i,2)=2*(x(i)^2)/(x(i)^2-(2*H)^2);
        term(i,3)=((x(i)-2*W)*x(i))/((x(i)-2*W)^2);
        term(i,4)=((x(i)+2*W)*x(i))/((x(i)+2*W)^2);
        term(i,5)=-2*(x(i)^2)/(x(i)^2-H^2);
        term(i,6)=-((x(i)+W)*x(i))/((x(i)+W)^2);
        term(i,7)=-((x(i)-W)*x(i))/((x(i)-W)^2);
        
        Gamma(i)=sum(term(i,:));
        
end

%now for the y part of the profile

y=x;


for i=1:size(x,1)
        
        term(i,1)=1;
        term(i,2)=((y(i)-2*H)*y(i))/((y(i)-2*H)^2);
        term(i,3)=((y(i)+2*H)*y(i))/((y(i)+2*H)^2);
        term(i,4)=2*(y(i)^2)/(y(i)^2-(2*W)^2);
        term(i,5)=-((y(i)-H)*y(i))/((y(i)-H)^2);
        term(i,6)=-((y(i)+H)*y(i))/((y(i)+H)^2);
        term(i,7)=-2*(y(i)^2)/(y(i)^2-W^2);
        
        Gamma(i)=sum(term(i,:));
end

figure(1); hold on; grid on; set(0,'DefaultAxesFontSize', 18)
h=plot(100*x/(W/2),Gamma,'b-');
axis([0 50 0 1.39]); title('Circulation vs percentage of tunnel width');
xlabel('100(r)/2W'); ylabel('Gamma(r)/GammaBW');
%legend('X direction','Y direction');
set(h(1),'LineWidth',2); grid minor;