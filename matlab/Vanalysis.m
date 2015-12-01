%==========================================================================
%           [Rcore,Haz,Hax]= Vanalysis(Run,strcore,strdyn,strmovie)
%========================================================================== 
% This function/script has been writen in order to process data sets for 
% axial vortex Particle Image Velocimetry data taken at the ODU lowspeed 
% windtunnel. It requires 3d vector data output by INSIGHT as well as
% averaged, pre-processed dynamic data output by the 'Load3d.m' script.
%
% INPUTS
% Run(#)     	  :Folder containing AvgData.xls file
% strcore(y/n)    :should be marked yes, or 'y' nearly always.
% strdyn(y/n)     :enables processing of Reynolds stress and Turbulence-E.
% strmovie(y/n)   :enables movie rendering of realtime data for given run,
%                   and takes a LONG time.
% Drive (char)    :drive letter for file paths, usually my 'E' drive.
%
%
% OUTPUTS 
% Rcore     	  :the calculated radius of the core
% Haz             :The maximum azimuthal velocity
% Hax             :The maximum axial velocity, usually slightly above Vfree
% [plots]         :Plots are all just saved into pre-specified directories, 
%                   the codes DOES NOT automatically create these 
%                   directories if they do not yet exist.
%
%========================================================================== 
%                           Function Syntax
%========================================================================== 

function [Rcore,Haz,Hax,Lax,Vaz,VazProfile,Vax,VaxProfile,Vra]=...
                               Vanalysis(Run,strcore,strdyn,strmovie,drive)
     close all;

%==========================================================================
%                  Secondary User defined parameters
%========================================================================== 
% File structure inputs, most to be used in script mode

    if exist('Run')==0      % If used as called function, then 'Run' will
    clear all; close all;   % already exist. If this m file is used as a 
    Run=60;                 % script, specify desired run number here.
    drive='E';           	% drive letter where data is being stored
    
    strcore='y';            % 'y' to manually enter core coordinates
    strdyn='y';             % 'y' to process Reynolds Stress and TurbEnergy
                            % must be 'y' for TEplots and RSplots to be on
    strmovie='n';           % 'y' to render real time dyndata movie. TIME++
    end    
     
    group=round(Run/10+0.4);% position grouping, 10 velocities per location

% general inputs
    rho=1.19;               % density of air
    r=25;                   % expected Max core radius in grid units. Used 
                            % to filter out edge noise when searching for 
                            % maximums near core
    
% data filtering inputs   
    Tol=0;                  % initial condition for Tol, variable to decide 
                            % how many adjacent spaces must be off by a 
                            % factor of (err) to flag a bad value. LEAVE 0.
    err=3;                  % If the Vp value of a point is this %/100 away 
                            % from at least Tol adjacent spaces, filter it 
    smoothcount =0;         % initial counter for number of smoothed vector
    repeat=1:3;             % number of smoothing iterations to run.
    
% Plot parameters
    TEplot= 'on';           % if on, turbulence energy plots will render
    RSplot= 'on';           % if on, Reynolds stress plots will render
    Vproff= 'on';           % if on, Velocity profile plot will render
    Vecfld= 'on';          % if on, vector field plot will render
    Vcont=  'on';          % if on,  V Contour plots will render
    
    SavePlots='y';          % if 'y', plots will be saved. OW
    SaveFigs='y';           % if 'y', figures will be saved. OW
    SaveProfData='y';       % if 'y', profile xls docs will be saved. OW
    
    bound= 70;              %(mm) bound on profile plots.
    
% output plot formats
    d=':\Data2\Oct_2014_Output\';  
                            %output filepath without drive letter
    form='-djpeg90';        %format of output (jpeg of 85% quality). type
                            %'help print' for more info
                            
    res='-r800';            %resolution to output charts, in ppi
    
% Define global varaibles
    global DynData;
    
%==========================================================================
%                 load Averaged data (pre-processed by 'Load3d.m')
%========================================================================== 

% open up all the processed data files   
    path=strcat(drive,':\Data2\Ely_May28th\Vector\',num2str(Run));cd(path);
    fid=fopen('AvgData.xls','r+');
    AvgData=xlsread('AvgData.xls');
    
% load Velocity index file
    fid=fopen(strcat(drive,':\Data2\Vindex.txt','r+'));
    Vin=importdata(strcat(drive,':\Data2\Vindex.txt'));
    Vfree=Vin(Run,2);   
    
% get rid of NaN Values, set them to zero (should no longer be needed)
  for a=1:size(AvgData,1)
     for b=1:size(AvgData,2)
         if isnan(AvgData(a,b))
            AvgData(a,b)=0;
         else end
     end
  end
    
%format of AvgData= [x,y,z,u,v,w,uU,vU,wU,N] where U=std and N=numpoints

%========================================================================== 
%       transforming tabulated data into 2 dimensional matrix of values 
%       and twin coordinate arrays
%========================================================================== 
% Transform Velocity data into a format matlab likes more.
disp('Transforming');

    k=1;j=1;
 for i=2:size(AvgData,1)
     if AvgData(i,2)==AvgData(i-1,2)
         j=j+1;
     else
         k=k+1; j=1;
         Y(k-1)=AvgData(i-1,2);
         Y(k)=AvgData(i,2);
     end
     
     %Vp = Velocity magnitudes in plane. (u,v components)
         Vp(k,j)=(AvgData(i,4).^2+AvgData(i,5).^2).^0.5;
     %VM= velocity magnitudes with all three components
         VM(k,j)=(AvgData(i,4).^2+AvgData(i,5).^2+AvgData(i,6)^.2).^0.5;
     %Velocities u,v,w in x,y,z direction.
         Vu(k,j)=AvgData(i,4);
         Vv(k,j)=AvgData(i,5);
         Vw(k,j)=AvgData(i,6);
     %Uncertainty in those velocities (Standard Deviation)
         Vuu(k,j)=AvgData(i,7);
         Vvu(k,j)=AvgData(i,8);
         Vwu(k,j)=AvgData(i,9);
 end
 
%==========================================================================
%    Smoothing out data: fills in small missing chunks or reduces spikes
%========================================================================== 
% check for bad values using simple difference approach, iterated (repeat) 
% times. when a sector is found which is (err) times greater in value 
% than at least (Tol) of the adjacent locations, it becomes an average of 
% all four locations around it. Large missing chunks will still cause
% issues.

for z=repeat
    for i=2:size(Vp,1)-1;
        for j=2:size(Vp,2)-1;
            if abs(Vp(i,j)-Vp(i,j+1))/Vp(i,j)>=err
                Tol=Tol+1; end 
            if abs(Vp(i,j)-Vp(i,j-1))/Vp(i,j)>=err
                Tol=Tol+1; end
            if abs(Vp(i,j)-Vp(i+1,j))/Vp(i,j)>=err
                Tol=Tol+1; end
            if abs(Vp(i,j)-Vp(i-1,j))/Vp(i,j)>=err
                Tol=Tol+1; end
            if Tol >=3
                Vp(i,j)=(Vp(i,j+1)+Vp(i,j-1)+Vp(i+1,j)+Vp(i-1,j))/4;
                VM(i,j)=(VM(i,j+1)+VM(i,j-1)+VM(i+1,j)+VM(i-1,j))/4;
                Vu(i,j)=(Vu(i,j+1)+Vu(i,j-1)+Vu(i+1,j)+Vu(i-1,j))/4;
                Vv(i,j)=(Vv(i,j+1)+Vv(i,j-1)+Vv(i+1,j)+Vv(i-1,j))/4; 
                Vw(i,j)=(Vw(i,j+1)+Vw(i,j-1)+Vw(i+1,j)+Vw(i-1,j))/4;
            %these values below are not correct uncertainty values,
            %but are adequate fillers until those values are needed.
                Vuu(i,j)=err; Vvu(i,j)=err; Vwu(i,j)=err;
                Tol=0;
                smoothcount=smoothcount+1;
            end
        end
    end
end

% alert user that vectors have been smoothed
    fprintf('%3.0f erant vectors have been removed\n',smoothcount);

% extract X and Y coordinates as arrays from the data matrix.
    X=AvgData(1:j,1); Y=Y';
    
% finding the maximum of Vp. Naming convention H=high (max), L=low (min).
    Hp=max(max(Vp'));
    [YHploc XHploc]=find(Vp==Hp);
    YHp=Y(YHploc); XHp=X(XHploc);
    
%========================================================================== 
%   Using in-plane velocities to estimate Vazimuthal in order to find core.
%========================================================================== 

% looking for a minimum wich is near the maximum, attempting to hone in on
% the vortex core on the first try. range (r) is an input value at top.

if exist('strcore')==0
    prompt='input manual X and Y index to find minimum? (y/n):';
    strcore = input(prompt,'s');
end
    if strcore=='n'
        %Works ONLY for the cleanest of data sets, manual option is best.
        [Lp]=min(min(Vp(YHploc-r:YHploc+r,XHploc-r:XHploc+r)'));
        [YLploc XLploc]=find(Vp==Lp);
        Yl=Y(YLploc);Xl=X(XLploc);
    else
        
        %If it exists, load manual core location from CoreIndex file.
        fopen(strcat(drive,':\Data2\CoreIndex.txt'))
        CoreIndex=importdata(strcat(drive,':\Data2\CoreIndex.txt'));
        
        if CoreIndex(Run,2)~=1
            XLploc=CoreIndex(Run,2);
            YLploc=CoreIndex(Run,3);
        else
            
                
        % output figures so the user may find the core themselves.
        tts('ALERT: User input required'); %text2speech
    
        figure(1); axis equal;
        quiver(AvgData(:,1),AvgData(:,2),AvgData(:,4),AvgData(:,5));
        xlabel('mm');ylabel('mm');title('Flowfield');grid on;
            Xin(:,2)=X(:,1);Xin(:,1)=1:length(X);
            prompt='Pease select from the list of index values of X:';
            XLploc=str2num(input(prompt,'s'));
            Yin(:,2)=Y(:,1);Yin(:,1)=1:length(Y);
            prompt='Pease select from the list of index values of Y:';
            YLploc=str2num(input(prompt,'s'));   
        end
        Lp=Vp(YLploc,XLploc);
    end
    
% getting a more precise core location, since the core need not lay on top 
% of a preexisting grid point. Examines the adjacent points to Min point
% Core Zone (cz)     CZ(:,:,#)===> 1=V, 2=X coord 3=Y coord
 
    CZ(:,:,1)=Vp(YLploc-2:YLploc+2,XLploc-2:XLploc+2);CZ(:,:,2:3)=zeros();
    for i=1:5
        for j=1:5
            CZ(i,j,2)=X(XLploc-3+j);
            CZ(i,j,3)=Y(YLploc-3+i);
        end
    end   

% if another value in CZ is lower than the center, issue a warning.   
    if min(min(CZ(:,:,1)))~=CZ(3,3,1);
        [aa,bb]=find(CZ==min(min(CZ(:,:,1))));
        XLploc=XLploc+bb-3;
        YLploc=YLploc+aa-3;
        Lp=Vp(YLploc,XLploc);
    fprintf('Core location amended to X= %3.0f Y= %3.0f \n',XLploc,YLploc);
    end
    
    % recalculate based on amended core location   
    CZ(:,:,1)=Vp(YLploc-2:YLploc+2,XLploc-2:XLploc+2);CZ(:,:,2:3)=zeros();
    for i=1:5
        for j=1:5
            CZ(i,j,2)=X(XLploc-3+j);
            CZ(i,j,3)=Y(YLploc-3+i);
        end
    end

% subtract the max value from velocity and take a weighted local average. 
% Allows us to find a more specific location of the core which almost
% certainly lies between grid ponts.
    A=max(max(CZ(:,:,1)));
    CZ(:,:,1)=CZ(:,:,1)-A;
    CZ(:,:,4)=(CZ(:,:,1).*-CZ(:,:,2));
    CZ(:,:,5)=(CZ(:,:,1).*-CZ(:,:,3));
    
    Xcore=-sum(sum(CZ(:,:,4)))/sum(sum(CZ(:,:,1)));
    Ycore=-sum(sum(CZ(:,:,5)))/sum(sum(CZ(:,:,1)));  

%==========================================================================
%                   finding real azimuthal velocities
%==========================================================================     
% Back to square one, now that the core has been found, we can calculate
% true and precise azimuthal velocities(Vaz) (ASUMING THE CORE IS CIRCULAR)
% this will need to be examined more closely if eliptical approximations
% prove necessary.

for xi=1:length(X)
    for yi=1:length(Y)
        %calculate azimuthal velocities (Vaz) and uncertainty (Vazu).
            A=[-(Y(yi)-Ycore);(X(xi)-Xcore);0]./...
                                (((Y(yi)-Ycore)^2+(X(xi)-Xcore)^2)^0.5);
            B=[Vu(yi,xi);Vv(yi,xi);0];
            Vaz(yi,xi)=abs((dot(A,B)));
            % Bu=[Vuu(yi,xi);Vvu(yi,xi);0];
            % Vazu(yi,xi)=abs((dot(A,Bu)));
        
        % calculate radial velocities (Vra) and uncertainty (Vrau).
            C=[(X(xi)-Xcore);(Y(yi)-Ycore);0]./...
                                (((Y(yi)-Ycore)^2+(X(xi)-Xcore)^2)^0.5);
            Vra(yi,xi)=-(dot(C,B));
            % Vrau(yi,xi)=abs((dot(C,Bu)));
        
        % transfer axial velocities(Vax) and uncertainty (Vaxu).
            Vax(yi,xi)=Vw(yi,xi);
            % Vaxu(yi,xi)=Vwu(yi,xi);
       
    end
end  

Lax=Vax(YLploc,XLploc);
    
% use Azimuthal velocity magnitudes to populate velocity profile function
ai=1; 
bi=1;
ci=1;
for xi=1:size(Vaz,2)
    for yi=1:size(Vaz,1)
        %only consider points which have a logical axial velocity,
        %indicating a likelyhood the vector is non-spurious
        if Vax(yi,xi)>=0.95*Vax(YLploc,XLploc)
            
        % only look at points near the core to help filter edge data
            if (((Y(yi)-Ycore)^2+(X(xi)-Xcore)^2)^0.5)<=bound
               
            %distance
            VazProfile(ai,1)=((Y(yi)-Ycore)^2+(X(xi)-Xcore)^2)^0.5;
            VaxProfile(ai,1)=VazProfile(ai,1);
                if X(xi)-Xcore<=0
                    VaxProfile(ai,1)=-VaxProfile(ai,1);
                end
                
            %complete velocity profile scatters
                VazProfile(ai,2)= Vaz(yi,xi);
                VaxProfile(ai,2)= Vax(yi,xi)/Vfree;
            
            %select vertical Y slice
                if abs(Xcore-X(xi))<=4
                    VaxProfileSlice(bi,:)=VaxProfile(ai,:);
                    bi=bi+1;
                end
                
                if abs(Xcore-X(xi))<=26
                    VazProfileSlice(ci,1:2)=VazProfile(ai,1:2);
                    ci=ci+1;
                end
            
                ai=ai+1;
            end
        end
    end
end
% organizes the values in ascending order and finds moving averages
    VazProfile=sortrows(VazProfile);
    VazProfileSlice=sortrows(VazProfileSlice);
    VazMovAvg(:,1) = VazProfileSlice(:,1);
    VazMovAvg(:,2) = smooth(VazProfileSlice(:,1),...
        VazProfileSlice(:,2),197,'rloess');
    VazMovAvg(:,2) = smooth(VazMovAvg(:,2),9);
    
    VaxProfile=sortrows(VaxProfile);
    VaxProfileSlice=sortrows(VaxProfileSlice);
    VaxMovAvg(:,1) = VaxProfileSlice(:,1);
    VaxMovAvg(:,2) = smooth(VaxProfileSlice(:,1),...
        VaxProfileSlice(:,2),29,'rloess');
    VaxMovAvg(:,2) = smooth(VaxMovAvg(:,2),9);

    xlswrite(strcat('VazProfile',num2str(Run)),VazMovAvg);
    
    xlswrite(strcat('VaxProfile',num2str(Run)),VaxMovAvg);
    
    %VaxProfindex=VaxMovAvg(:,1);
    %VaxProf=(VaxMovAvg(:,2));
    
% create Vtheta*r data for some circulation indication
    CProfile=VazProfile;
    CProfile(:,2)=CProfile(:,1).*CProfile(:,2)*2*pi()/1000;
    
% inverting the profile (IP) to find the few max velocities.
    VazIP(:,1)=VazProfile(:,2); VazIP(:,2)=VazProfile(:,1); 
    VazIP=flipud(sortrows(VazIP));
    
% IMPORTANT: this section has accuracy implications!!! 
% Rather than relying on the one single maximum value to indicate the core 
% location, it is wise to consider other points which are near the maximum
% to yield more averaged and arguably accurate core size readings. the
% number of top values to consider is up for debate, 5-20 works well.
        Haz=mean(VazIP(1:5,1));
        Rcore=mean(VazIP(1:5,2));

% finding the maximum of Vradial
    Hra=max(max(Vra'));
    [YHraloc XHraloc]=find(Vra==Hra);
    YHra=Y(YHraloc); XHra=X(XHraloc);
%finding the minimum of Vradial
    Lra=min(min(Vra'));
    [YLraloc XLraloc]=find(Vra==Lra);
    YLra=Y(YLraloc); XLra=X(XLraloc);
%finding the maximum of Vaxial (Vax)
    Hax=max(max(Vax'));
    [YHaxloc XHaxloc]=find(Vax==Hax);
    YHax=Y(YHaxloc); XHax=X(XHaxloc); 
    
%========================================================================== 
%               theoretical plot data preparation
%==========================================================================  
   
%Ash Zuckerwar / Scully Kaufmann vortex
   %Vtheta(r)=2*Vtheta_max *(r/Rcore)/((r/Rcore)^2+1)
   rt=0:.05:bound;
   VazAZ(:,1)=rt;
   VazAZ(:,2)=2*Haz*(rt./Rcore)./(((rt./Rcore).^2)+1);
    
%rankine vortex
   %Vtheta(r)=gamma/2*pi*r for r>R
   %Vtheta(r)=gamma*r/(2*pi*R^2) r<R
   
   VazRankine(:,1)=rt;
    parfor i=1:length(rt)
        if rt(i)<= Rcore
            VazRankine(i,2)=Haz*rt(i)/Rcore;
        else
            VazRankine(i,2)=Haz*Rcore/rt(i);
        end
    end
    
    %calculate core circle for plotting purposes.
    ang=0:0.01:2*pi; 
    CircX=Rcore*cos(ang);
    CircY=Rcore*sin(ang);
    CircZ(1:length(CircX))=Haz;
    
    %Core radius line equation for plotting purposes
    RcoreMarker(1,:)=[Rcore 0]; 
    RcoreMarker(2,:)=[Rcore max(CProfile(:,2))];
 
%========================================================================== 
%                       creating output data
%==========================================================================     
      
    disp('====================================================');
    fprintf('Run # %2.1f \nFree Stream Velocity = %2.2f m/s\n',Run,Vfree);
    disp('====================================================');
    fprintf('Max axial velocity found to be %2.2f m/s\n',Hax);
    fprintf('Core found at \t X= %3.4f, \tY= %3.4f \n',Xcore,Ycore);
    fprintf('Core raduis = %3.4f (mm) \n',Rcore);
    fprintf('Vtheta_max = %3.4f (m/s) \n',Haz);
    
%Make pretty pictures
    
%plot AvgData as a vector field plot
if Vecfld(1:2)=='on'
figure(1); axis equal;grid on;
    t= strcat('Flowfield Plot [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');
    quiver(AvgData(:,1),AvgData(:,2),AvgData(:,4),AvgData(:,5));
    xlabel('mm');ylabel('mm');title(t);
end

if Vproff(1:2)=='on'
   
%Circulation profile
figure (6); hold on; grid on;
    t= strcat('Circulation scatter [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    h=plot(CProfile(:,1),CProfile(:,2),'b.',...
        RcoreMarker(:,1),RcoreMarker(:,2),'r-');
    set(h(1),'markersize',3);
    set(h(2),'LineWidth',1);
    legend('discrete point',strcat('Core Radius = ',num2str(Rcore),'mm'));
    legend('Location','SouthEast');
    xlabel('distance to core (mm)'); ylabel('Circulation Strength ');
    set(gca,'layer','top');
    
% Azimuthal velocity profile

        %change coordinates of core marker
            RcoreMarker(1,2)=0; 
            RcoreMarker(2,2)=Haz*1.15;
            
figure (5); hold on; grid on;
    t= strcat('Azimuthal velocity (v_a_z) scatter [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    axis([0 bound 0 Haz*1.15])
    h=plot(VazProfile(:,1),VazProfile(:,2),'b.',...
        RcoreMarker(:,1),RcoreMarker(:,2),'r-'); 
        %,VazMovAvg(:,1),VazMovAvg(:,2),'k-');
    set(h(1),'markersize',3);
    set(h(2),'LineWidth',1);
    %set(h(3),'LineWidth',1.5);
    legend('Location','NorthEast');
    legend('discrete point',strcat('Core Radius = ',num2str(Rcore),'mm'))
        %,...'Moving Average');
    xlabel('distance to core (mm)'); ylabel('v_a_z (m/s)');
    set(gca,'layer','top'); 
    
% Axial velocity profile   

        %change coordinates of core marker so it fits on new plot
            RcoreMarker(1,2)=0; 
            RcoreMarker(2,2)=1.1;
    
figure (7); hold on; grid on;
    t= strcat('Axial velocity (v_z) profile and scatter [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    h=plot(VaxProfile(:,1),VaxProfile(:,2),'b.',...
        RcoreMarker(:,1),RcoreMarker(:,2),'r-',...
        -RcoreMarker(:,1),RcoreMarker(:,2),'r-');
        %,...VaxMovAvg(:,1),VaxMovAvg(:,2),'k-');  
    set(h(1),'markersize',3);
    set(h(2:3),'LineWidth',1);
    %set(h(4),'LineWidth',1.5);
    axis([-bound bound 0.7 1.1]);
    legend('Vax points',strcat('Core Radius = ',num2str(Rcore),'mm'));
    legend('Location','SouthEast');
    xlabel('distance to core (mm)'); ylabel('v_z (%V_free)');
    set(gca,'layer','top');    
end   
   

if Vcont(1:2)=='on'

%Azimuthal velocity contour map
figure (2); hold on;axis equal; 
    t= strcat('Azimuthal velocity (v_az) [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(Vaz,size(Vaz,1)*size(Vaz,2),1); top=top(top~=-0.01);
    top=mean(top)+4*std(top);
    scale=[0:Haz/20:Haz];
    [cs,h]=contourm(Y,X,Vaz,scale); 
    leg=clegendm(cs,h,-1,' m/s');
    contourf(X,Y,Vaz,scale); contour(X,Y,Vaz,scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')
    
%Radial velocity contour map    
figure (3); hold on;axis equal; colormap(RedBlue);
    t= strcat('Radial velocity (v_r) (-) is outwards [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(Vra,size(Vra,1)*size(Vra,2),1); top=top(top~=-0.01);
    top=mean(top)+4*std(top);
    scale=[-top:top/8:top];
    [cs,h]=contourm(Y,X,Vra,scale); 
    leg=clegendm(cs,h,-1,' m/s');
    contourf(X,Y,Vra,scale); contour(X,Y,Vra,scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')
    
%Axial velocity contour map    
figure (4); hold on;axis equal; colormap(jet);
    t= strcat('Axial velocity (v_z) [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    Low=Vax(YLploc,XLploc)*.9;
    scale=[Low:(Vfree*1.2-Low)/20:Vfree*1.2];
    [cs,h]=contourm(Y,X,Vax,scale); 
    leg=clegendm(cs,h,-1,' m/s');
    contourf(X,Y,Vax,scale); contour(X,Y,Vax,scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top') 

end

%========================================================================== 
%           dynamic display and non-averaged value calculations
%========================================================================== 
%dynamic display (individual plots of all 200 runs if desired)

if exist('strdyn')==0
    prompt='begin processing for dynamic data and turb intensity? (y/n):';
    strdyn = input(prompt,'s');
end
    
if strdyn=='y'
    disp('Loading');
    H=csvread('DynamicData.txt');
    disp('Transforming');
    DynData=reshape(H,length(H),6,size(H,2)/6);
    
%Transform Dynamic velocity data in the same way as before, put them in
%terms of azimuthal and radial velocities.

for m=1:size(DynData,3)
    
    k=1;j=1;
    for i=2:size(DynData,1)
        if DynData(i,2,m)==DynData(i-1,2,m)
            j=j+1;
        else
            k=k+1; j=1;
        end
        
        if j <= size(X,1)
        %create additional data maps 
            DynVu(k,j,m)=DynData(i,4,m);
            DynVv(k,j,m)=DynData(i,5,m);
            DynVw(k,j,m)=DynData(i,6,m);
        %since we already know the core location, lets put in cylindrical
        %azimuthal
            A=[-(Y(k)-Ycore);(X(j)-Xcore);0]./...
                                (((Y(k)-Ycore)^2+(X(j)-Xcore)^2)^0.5);
            B=[DynVu(k,j,m);DynVv(k,j,m);0];
            DynVaz(k,j,m)=abs((dot(A,B)));
        
        %Radial
            C=[(X(j)-Xcore);(Y(k)-Ycore);0]./...
                                (((Y(k)-Ycore)^2+(X(j)-Xcore)^2)^0.5);
            DynVra(k,j,m)=-(dot(C,B));
        end
    end
end

disp('Calculating Turbulence and Reynolds stress Values');
for k=1:size(DynVu,1)
    for j=1:size(DynVu,2)
        
%reynolds stresses in the 3 directional combinations
        Ar= abs(DynVra(k,j,:)-Vra(k,j)); 
        Br=(Ar(Ar~=0)).^2; 
        Cr=sqrt(sum(Br)/length(Br));
        
        At= abs(DynVaz(k,j,:)-Vaz(k,j)); 
        Bt=(At(At~=0)).^2; 
        Ct=sqrt(sum(Bt)/length(Bt));
        
        Aw= abs(DynVw(k,j,:)-Vw(k,j)); 
        Bw=(Aw(Aw~=0)).^2; 
        Cw=sqrt(sum(Bw)/length(Bw));
        
        RSrt(k,j)=Cr*Ct/(Vfree*Vfree); 
            if isnan(RSrt(k,j))==1
                RSrt(k,j)=0;
            end
        RSrw(k,j)=Cr*Cw/(Vfree*Vfree);
            if isnan(RSrw(k,j))==1
                RSrw(k,j)=0;
            end
        RStw(k,j)=Ct*Cw/(Vfree*Vfree);
            if isnan(RStw(k,j))==1
                RStw(k,j)=0;
            end
        
% turbulent energies (changed to reflect reynolds stress format) and stored
% in one 3 dimensional matrix

        if VM(k,j)>=.05
            TI(k,j,2)= (Cr^2)/(Vfree*Vfree);  %radial
            TI(k,j,1)= (Ct^2)/(Vfree*Vfree);  %azimuthal
            TI(k,j,3)= (Cw^2)/(Vfree*Vfree);  %axial
        else 
            TI(k,j,1:3)=0;
        end
    end
end



%========================================================================== 
%                        Core Wandering analysis
%========================================================================== 
% examine the variance tied to tangential velocity in the core area to
% account for wandering as described by Davenport (1996). This analysis
% assumes that ALL the variance in azimuthal velocity at the core location
% is rooted in core wandering, so this data means nothing without a general
% estimate of azimuthal velocity errors.

    % XLploc
    % YLploc
    
    Wand = TI(YLploc,XLploc,1);
    Wand = (Wand*(Vfree*Vfree))^.5;
    TrueMeasure =(1-(2*1.25643*Wand*Wand/((Rcore)^2)))^0.5;
    fprintf('True radius / measured radius = %1.3f', TrueMeasure);
    
%========================================================================== 
%       Real time Movie Generation for VERY complete data sets
%========================================================================== 

if exist('strmovie')==0
    prompt='Render video of velocity fluctuations? (y/n):';
    strmovie = input(prompt,'s');
end

if strmovie=='y'

%for the sake of easier viewing and continuity, setting absent values to
%zero makes it tough to distinguish between a vortex core and an absent
%value. So we set 0 values to -10, this seems to work well visually.

    for i=1:size(DynVaz,1)
        for j=1:size(DynVaz,2)
            for m=1:size(DynVaz,3)
                if DynVaz(i,j,m)==0
                    DynVazDiff(i,j,m)=-10;
                else end
            end
        end
    end
    
%Subtract the average in plane velocities to show deviations from avg.
    
for i=1:size(DynVaz,1)
    for j=1:size(DynVaz,2)
        for k=1:size(DynVaz,3)
            if DynVazDiff(i,j,k)~=-10
                DynVazDiff(i,j,k)=DynVaz(i,j,k)-Vaz(i,j);
            end
        end
    end
end
                

%jumpstart it by oppening the first figure outside the loop, for some
%reason matlab has a nervous breakdown if you don't do this.
    
figure(8);axis equal; 
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    pause(1);
    scale=[-Haz/6:Haz/24:Haz/6];
    
%the rest of the frames, after jumpstarting the first frame

for i=1:size(DynVaz,3)
    title('Azimuthal velocity deviation from average');
    contourf(X,Y,DynVazDiff(:,:,i),scale); hold on;
    contour(X,Y,DynVazDiff(:,:,i),scale);
    
    [cs,h]=contourm(Y,X,DynVazDiff(:,:,i),scale); 
    leg=clegendm(cs,h,1,'m/s');
    plot3(Xcore,Ycore,Lp,'og');
    plot(Xcore+CircX,Ycore+CircY,'-g');
    grid on; set(gca,'layer','top'); axis equal; grid on;
    Movie(i)=getframe; hold off;

end

    %write AVI file to HDD, these are usually quite large
    disp('writing animation to .avi file');
    movie2avi(Movie,'IP_Velocitydiff3fps','FPS',2.87);
end

%==========================================================================   
%     Plots plots plots plots plots plots plots plots plots plots plots 
%     eeeeerybody!....
%========================================================================== 

% Turbulence energy maps: (may be a misnomer at this point) 
% note that these are component wise (v")^2 / Vfree^2 plots.

if TEplot(1:2)=='on'
figure (9); hold on;axis equal;
    t= strcat('RS Azimuthal v_a_z^{"2}/U^{2} [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(TI(:,:,1),size(TI,1)*size(TI,2),1); top=top(top~=0);
    top=mean(top)+5*std(top);
    scale=[0:top/20:top];
    [cs,h]=contourm(Y,X,TI(:,:,1),scale); 
    leg=clegendm(cs,h,-1);
    contourf(X,Y,TI(:,:,1),scale); contour(X,Y,TI(:,:,1),scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')

figure (10); hold on;axis equal;
    t= strcat('RS Azimuthal v_r^{"2}/U^{2} [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(TI(:,:,2),size(TI,1)*size(TI,2),1); top=top(top~=0);
    top=mean(top)+1*std(top);
    scale=[0:top/20:top];
    [cs,h]=contourm(Y,X,TI(:,:,2),scale); 
    leg=clegendm(cs,h,-1);
    contourf(X,Y,TI(:,:,2),scale); contour(X,Y,TI(:,:,2),scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')

figure (11); hold on;axis equal;
    t= strcat('RS Azimuthal v_z^{"2}/U^{2} [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(TI(:,:,3),size(TI,1)*size(TI,2),1); top=top(top~=0);
    top=mean(top)+1*std(top);
    scale=[0:top/20:top];
    [cs,h]=contourm(Y,X,TI(:,:,3),scale); 
    leg=clegendm(cs,h,-1);
    contourf(X,Y,TI(:,:,3),scale); contour(X,Y,TI(:,:,3),scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')
end

%  Reynolds Stresses
if RSplot(1:2)=='on'
figure (12); hold on;axis equal;
    t= strcat('RS v_r^{"}*v_z^{"}/U^{2} [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(RSrw,size(RSrw,1)*size(RSrw,2),1); top=top(top~=0);
    top=mean(top)+3*std(top);
    scale=[0:top/20:top];
    [cs,h]=contourm(Y,X,RSrw,scale); 
    leg=clegendm(cs,h,-1);
    contourf(X,Y,RSrw,scale); contour(X,Y,RSrw,scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')

figure (13); hold on;axis equal;
    t= strcat('RS v_r^{"}*v_a_z^{"}/U^{2}  [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(RSrt,size(RSrt,1)*size(RSrt,2),1); top=top(top~=0);
    top=mean(top)+3*std(top);
    scale=[0:top/20:top];
    [cs,h]=contourm(Y,X,RSrt,scale); 
    leg=clegendm(cs,h,-1);
    contourf(X,Y,RSrt,scale); contour(X,Y,RSrt,scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')

figure (14); hold on;axis equal;
    t= strcat('RS v_z^{"}*v_a_z^{"}/U^{2}  [Run #',...
        num2str(Run),' : Position',...
        num2str(group),' : U= ',num2str(Vfree),'m/s ]');title(t);
    top=reshape(RStw,size(TI,1)*size(RStw,2),1); top=top(top~=0);
    top=mean(top)+3*std(top);
    scale=[0:top/20:top];
    [cs,h]=contourm(Y,X,RStw,scale); 
    leg=clegendm(cs,h,-1);
    contourf(X,Y,RStw,scale); contour(X,Y,RStw,scale);
    plot3(Xcore,Ycore,Lp,'ow');
    plot(Xcore+CircX,Ycore+CircY,'-w');
    grid on; set(gca,'layer','top')
end

end

%==========================================================================
%   saving all of those plots as jpeg images to be included in thesis
%========================================================================== 

if SaveProfData=='y'
    disp('Exporting Velocity profiles xls tables');
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    xlswrite(strcat('VazProfile',num2str(Run),'.xls'),VazMovAvg);
    xlswrite(strcat('VaxProfile',num2str(Run),'.xls'),VaxMovAvg);
end


if SavePlots=='y'
disp('Saving plots and charts!');

%Velocity plots
if Vcont(1:2)=='on'
   path=strcat(drive,d,'Velocity Azimuthal\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (2),strcat('Azimuthal V',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (2),strcat('Azimuthal V',num2str(Run)),form,res);
    
   path=strcat(drive,d,'Velocity Radial\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (3),strcat('Radial V',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (3),strcat('Radial V',num2str(Run)),form,res);
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
   path=strcat(drive,d,'Velocity Axial\',num2str(group));
   
    print(figure (4),strcat('Axial V',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (4),strcat('Axial V',num2str(Run)),form,res);

    
end   

%Profile plots, velocities and circulations etc.
if Vproff(1:2)=='on'
   path=strcat(drive,d,'Profile Azimuthal Velocity\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (5),strcat('Vaz Profile',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (5),strcat('Vaz Profile',num2str(Run)),form,res);
    
   path=strcat(drive,d,'Profile Circulation\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (6),strcat('Circulation Profile',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (6),strcat('Circulation Profile',num2str(Run)),form,res);
    
   path=strcat(drive,d,'Profile Axial Velocity\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (7),strcat('Vax Profile',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (7),strcat('Vax Profile',num2str(Run)),form,res);

end

if strdyn=='y'
%Turbulent energies
if TEplot(1:2)=='on'
    
   path=strcat(drive,d,'Turbulent Energy Azimuthal\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (9),strcat('Azimuthal TE',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (9),strcat('Azimuthal TE',num2str(Run)),form,res);

   path=strcat(drive,d,'Turbulent Energy Radial\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (10),strcat('Radial TE',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (10),strcat('Radial TE',num2str(Run)),form,res);

   path=strcat(drive,d,'Turbulent Energy Axial\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (11),strcat('Axial TE',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (11),strcat('Axial TE',num2str(Run)),form,res);
end

%Reynolds Stresses
if RSplot(1:2)=='on'
            
   path=strcat(drive,d,'Reynolds Stress Radial-Axial\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (12),strcat('Radial-Axial RS',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (12),strcat('Radial-Axial RS',num2str(Run)),form,res);

   path=strcat(drive,d,'Reynolds Stress Radial-Azimuthal\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (13),strcat('Radial-Aziumthal RS',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (13),strcat('Radial-Aziumthal RS',num2str(Run)),form,res); 

   path=strcat(drive,d,'Reynolds Stress Axial-Azimuthal\',num2str(group));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (14),strcat('Axial-Aziumthal RS',num2str(Run)),form,res);
    path=strcat(drive,d,num2str(Run));
    if exist(path)~=7
        mkdir(path);
    end
    cd(path);
    print(figure (14),strcat('Axial-Aziumthal RS',num2str(Run)),form,res);
end
end
end

%========================================================================== 
%                              Save figures
%========================================================================== 

if SaveFigs=='y'
disp('Saving Figures!');
path=strcat(drive,d,'Fig Files\');
if exist(path)~=7
    mkdir(path);
end
cd(path);

%Velocity plots
if Vcont(1:2)=='on'

    saveas(figure(2),strcat('Azimuthal V',num2str(Run)), 'fig');
    saveas(figure(3),strcat('Radial V',num2str(Run)), 'fig');
    saveas(figure(4),strcat('Axial V',num2str(Run)), 'fig'); 

end   
if Vproff(1:2)=='on'

    saveas(figure(5),strcat('Vaz Profile',num2str(Run)), 'fig'); 
    saveas(figure(6),strcat('Circulation Profile',num2str(Run)), 'fig'); 
    saveas(figure(7),strcat('Vax Profile',num2str(Run)), 'fig'); 
end

if strdyn=='y'
    
%Turbulent energies

if TEplot(1:2)=='on'
    
    saveas(figure(9),strcat('Azimuthal TE',num2str(Run)), 'fig'); 
    saveas(figure(10),strcat('Radial TE',num2str(Run)), 'fig'); 
    saveas(figure(11),strcat('Axial TE',num2str(Run)), 'fig'); 
end

%Reynolds Stresses

if RSplot(1:2)=='on'
            
    saveas(figure(12),strcat('Radial-Axial RS',num2str(Run)), 'fig'); 
    saveas(figure(13),strcat('Radial-Azimuthal RS',num2str(Run)), 'fig'); 
    saveas(figure(14),strcat('Axial-Azimuthal RS',num2str(Run)), 'fig'); 
end
end
end
%========================================================================== 
%                              Clean up
%========================================================================== 

% Reset default path from the beginning of the function.
path='E:\Data2'; cd(path);

disp(strcat('Run ',num2str(Run),' Complete, Hooorayyyyy!'));

%     strcore          
%     strdyn                   
%     strmovie       



end













