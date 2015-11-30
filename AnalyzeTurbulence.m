% Loads up the workspace containing Data for turbulent analysis.
clc; clear all; %close all;
addpath('C:/Users/Jeff/Desktop/Github/thesis-pivpr')
%Plot tuning parameters.
    Runs=[1:38 41:70];  %Runs to plot.
    Chord = 101.6;      %Chord length of model (in milimeters) (4 inches)
    n=1024;             %number of distinct colors on colormap
    c=20;               %number of distinct colors to display on each plot.  
    Lloc='NorthEast';   %location on chart to place maximum value.
    
    Xloc=[21.5 27.875 31 34 36 38 40].*(25.4)./Chord; %(milimeters)
    Xloc=round(10*Xloc)./10;
    
    
% Loading required data
     %load('CartesianDynamic.mat');
     load('CartesianRS.mat');
     load('CartesianTE.mat');

     %load('CylindricalDynamic.mat');
     load('CylindricalRS.mat');
     load('CylindricalTE.mat');

     load('Coordinates.mat'); 
     load('Rcore.mat'); 
     load('Vfree.mat');
     %load('DynamicVelocities.mat');
     
     cd('C:\Users\Jeff\Desktop\Github\thesis-pivpr\plots\2015_11_30');
 
for i=Runs
    close all;
    X{i}=X{i}-Xcore(i); X{i}=X{i}/Chord;
    Y{i}=Y{i}-Ycore(i); Y{i}=Y{i}/Chord;
    area=([min(X{1}) max(X{i}) min(Y{i}) max(Y{i})]);
    
    uuBar=abs(uuBar);
    vvBar=abs(vvBar);
    wwBar=abs(wwBar);
    uvBar=abs(uvBar);
    uwBar=abs(uwBar);
    vwBar=abs(vwBar);
    
    rrBar=abs(rrBar);
    ttBar=abs(ttBar);
    rtBar=abs(rtBar);
    rwBar=abs(rwBar);
    twBar=abs(twBar);
    
%--------------------------------------------------------------------------
%                           Cartesian coordinates
%--------------------------------------------------------------------------

    %find paddings to avoid edge noise contamination of maximums.
    lp=min(find(round(10*X{i})./10==-0.7));   %left padding
    rp=max(find(round(10*X{i})./10== 0.7));   %right padding

    top1=max(max(uvBar(:,lp:rp,i)));
    
    top2=max([max(max(uuBar(:,lp:rp,i))),max(max(vvBar(:,lp:rp,i))),...
        max(max(uvBar(:,lp:rp,i))),max(max(uwBar(:,lp:rp,i))),...
        max(max(vwBar(:,lp:rp,i)))]);
    
    top3= max(max(wwBar(:,lp:rp,i)));
    
    %set up the figure and manipulate the colormap
    hfig=figure (i); set(hfig,'units','normalized','outerposition',[0 0 1 1])
    a=suptitle(['Time Averaged Reynolds Stresses (Cartesian), Case',num2str(i),' [ z/c = ',...
        num2str(Xloc(round(i/10+0.49))),' , V_f_r_e_e = ',num2str(Vfree(i)),'m/s ]']);
    set(a,'Fontsize',12);
    subplot = @(m,n,p) subtightplot(m,n,p,[0 0.04], [0.06 0.04], [0.05 0.10]);
    
    col=colormap(cbrewer('seq','YlOrRd',n-round((n*top2/top3))));
    col=[flipud(colormap(cbrewer('seq','YlGnBu',round(0.5+(n*(top2-top1)/top3)))));col];
    col=[colormap(cbrewer('seq','Blues',round(0.5+(n*top1/top3))));col];
    colormap(col);
    
    set(get(gca,'xlabel'),'fontsize',6)
    
    ax(1)=subplot(2,3,3);hold on;
    [CS,H]=contourf(X{i},Y{i},wwBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none'); xlabel('x/c');
    h_title=title('(w''w'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area); xlabel('x/c');
    b=legend(strcat('Max = ',num2str(max(max(wwBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    set(get(gca,'xlabel'),'fontsize',6);
    xlabel('x/c');
    
    ax(2)=subplot(2,3,1);hold on;
    [cs,h]=contourf(X{i},Y{i},uuBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none'); 
    h_title=title('(u''u'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area); xlabel('x/c');
    b=legend(strcat('Max = ',num2str(max(max(uuBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');ylabel('y/c');
    
    ax(3)=subplot(2,3,2);hold on;
    [cs,h]=contourf(X{i},Y{i},vvBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(v''v'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area); xlabel('x/c');
    b=legend(strcat('Max = ',num2str(max(max(vvBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');

    ax(4)=subplot(2,3,4);hold on;
    [cs,h]=contourf(X{i},Y{i},uvBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(u''v'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(uvBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');ylabel('y/c');
    
    ax(5)=subplot(2,3,5);hold on;
    [cs,h]=contourf(X{i},Y{i},uwBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(u''w'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(uwBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');
    
    ax(6)=subplot(2,3,6);hold on;
    [cs,h]=contourf(X{i},Y{i},vwBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(v''w'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(vwBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6);legend(gca,'boxoff'); 
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');
    
    
    
    %moves the colorbar to the right hand side.
        h=colorbar;
        set(h,'Position',[.92 .1 .02 0.80]);

        for ii=[1 2 3 4 5 6]
            pos=get(ax(ii),'Position');
            set(ax(ii),'Position',[pos(1) pos(2) pos(3) pos(4)]);
            ticks=[0,top1,(top1+((top2-top1)/4)):((top2-top1)/4):top2,...
                (top2+((top3-top2)/12)):((top3-top2)/12):(top3-((top3-top2)/12)),top3*.999];
            ticks=round(1000*ticks)./1000;
            set(h,'YTick',ticks,'FontSize',5);
        end
        hold off;
        
    %Save the figure
        print(figure (i),strcat('RS_Run_',num2str(i),' Cartesian'),'-djpeg85','-r600');

 
%-------------------------------------------------------------------------- 
%                           Cylindrical coordinates
%--------------------------------------------------------------------------
    close all;
    
    top1=max(max(rtBar(:,lp:rp,i)));
    
    top2=max([max(max(rrBar(:,lp:rp,i))),max(max(ttBar(:,lp:rp,i))),...
        max(max(rwBar(:,lp:rp,i))),max(max(twBar(:,lp:rp,i)))]);
    
    %set up the figure and manipulate the colormap
    hfig=figure (i); set(hfig,'units','normalized','outerposition',[0 0 1 1])
    a=suptitle(['Time Averaged Reynolds Stresses (Cylindrical), Case',num2str(i),' [ z/c = ',...
        num2str(Xloc(round(i/10+0.49))),' , V_f_r_e_e = ',num2str(Vfree(i)),'m/s ]']);
    set(a,'Fontsize',12);
    subplot = @(m,n,p) subtightplot(m,n,p,[0 0.04], [0.06 0.04], [0.05 0.10]);
    
    col=colormap(cbrewer('seq','YlOrRd',n-round((n*top2/top3))));
    col=[flipud(colormap(cbrewer('seq','YlGnBu',round(0.5+(n*(top2-top1)/top3)))));col];
    col=[colormap(cbrewer('seq','Blues',round(0.5+(n*top1/top3))));col];
    colormap(col);
    
    set(get(gca,'xlabel'),'fontsize',6)
    
    ax(1)=subplot(2,3,3);hold on;
    [CS,H]=contourf(X{i},Y{i},wwBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(w''w'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(wwBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');
    
    ax(2)=subplot(2,3,1);hold on;
    [cs,h]=contourf(X{i},Y{i},rrBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(r''r'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(rrBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');ylabel('y/c');
    
    ax(3)=subplot(2,3,2);hold on;
    [cs,h]=contourf(X{i},Y{i},ttBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(t''t'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area); 
    b=legend(strcat('Max = ',num2str(max(max(ttBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');
    ax(4)=subplot(2,3,4);hold on;
    [cs,h]=contourf(X{i},Y{i},rtBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(r''t'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(rtBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');ylabel('y/c');
    
    ax(5)=subplot(2,3,5);hold on;
    [cs,h]=contourf(X{i},Y{i},rwBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(r''w'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(rwBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6); legend(gca,'boxoff');
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');
    
    ax(6)=subplot(2,3,6);hold on;
    [cs,h]=contourf(X{i},Y{i},twBar(1:length(Y{i}),1:length(X{i}),i),c,...
        'LineColor','none');
    h_title=title('(t''w'')'); set(h_title,'FontSize',10);
    grid on; set(gca,'layer','top'); axis(area);
    b=legend(strcat('Max = ',num2str(max(max(twBar(:,:,i))))));
    a=legend('Location',Lloc); set(a,'FontSize',6);legend(gca,'boxoff'); 
    caxis manual; caxis([0 top3]); set(gca,'FontSize',6); axis equal
    xlabel('x/c');
    
    
    
    %moves the colorbar to the right hand side.
        h=colorbar;
        set(h,'Position',[.92 .1 .02 0.80]);

        for ii=[1 2 3 4 5 6]
            pos=get(ax(ii),'Position');
            set(ax(ii),'Position',[pos(1) pos(2) pos(3) pos(4)]);
            ticks=[0,top1,(top1+((top2-top1)/4)):((top2-top1)/4):top2,...
                (top2+((top3-top2)/12)):((top3-top2)/12):(top3-((top3-top2)/12)),top3*.999];
            ticks=round(1000*ticks)./1000;
            set(h,'YTick',ticks,'FontSize',5);
        end
        hold off;
        
    %Save the figure

        print(figure (i),strcat('RS_Run_',num2str(i),' Cylindrical'),'-djpeg85','-r600');
end


% this section was used to atually create the giant matrixes of data

% BatchCropDirectory('E:\Data2\Ely_May28th\New Output',3,'jpg');

% for i =1:70
%     t{i}=zeros(length(Y{i}),length(X{i}),200);
%     r{i}=zeros(length(Y{i}),length(X{i}),200);
%     fprintf('Run %2.0f zero filled \n',i);
%     
%     for xi = 1:length(X{i})
%         fprintf('At %2.3f%% of Total completion \n',(100/70)*(xi/length(X{i}))+100*(i-1)/70);
%         for yi = 1:length(Y{i})
%             for l =1:size(u{i},3)
%              
%                 A=[-(Y{i}(yi)-Ycore(i));(X{i}(xi)-Xcore(i));0]./...
%                     (((Y{i}(yi)-Ycore(i))^2+(X{i}(xi)-Xcore(i))^2)^0.5);
%                 C=[(X{i}(xi)-Xcore(i));(Y{i}(yi)-Ycore(i));0]./...
%                     (((Y{i}(yi)-Ycore(i))^2+(X{i}(xi)-Xcore(i))^2)^0.5);
%                 B=[u{i}(yi,xi,l);v{i}(yi,xi,l);0];
% 
%                 t{i}(yi,xi,l)=(dot(A,B));   %Azimuthal
%                 r{i}(yi,xi,l)=-(dot(C,B));  %Radial   
%             end
%         end
%     end
%     fprintf('Run %2.0f done \n',i);
% end

% for i =1:70
%     for j = 1:size(u{i},1)
%         for k = 1:size(u{i},2)
%             for l=1:size(u{i},3)
%                rr{i}(j,k,l)=r{i}(j,k,l)^2; 
%             end
%             
%             rrBar(j,k,i)=mean(rr{i}(j,k,:));
%         end
%     end
%     fprintf('Run %2.0f of rr done \n',i); 
% end
% 
% clear rr
% disp('Finished with rr'); 
% 
% for i =1:70
%     for j = 1:size(u{i},1)
%         for k = 1:size(u{i},2)
%             for l=1:size(u{i},3)
%                tt{i}(j,k,l)=t{i}(j,k,l)^2; 
%             end
%             
%             ttBar(j,k,i)=mean(tt{i}(j,k,:));
%         end
%     end
%     fprintf('Run %2.0f of tt done \n',i); 
% end
% 
% clear tt
% disp('Finished with tt'); 
% 
% for i =1:70
%     for j = 1:size(u{i},1)
%         for k = 1:size(u{i},2)
%             for l=1:size(u{i},3)
%                ww{i}(j,k,l)=w{i}(j,k,l)^2; 
%             end
%             
%             wwBar(j,k,i)=mean(ww{i}(j,k,:));
%         end
%     end
%     fprintf('Run %2.0f of ww done \n',i); 
% end
% 
% clear ww
% disp('Finished with ww'); 
% 
% for i =1:70
%     for j = 1:size(u{i},1)
%         for k = 1:size(u{i},2)
%             for l=1:size(u{i},3)
%                rt{i}(j,k,l)=r{i}(j,k,l)*t{i}(j,k,l); 
%             end
%             
%             rtBar(j,k,i)=mean(rt{i}(j,k,:));
%         end
%     end
%     fprintf('Run %2.0f of rt done \n',i); 
% end
% 
% clear rt
% disp('Finished with rt'); 
% 
% for i =1:70
%     for j = 1:size(u{i},1)
%         for k = 1:size(u{i},2)
%             for l=1:size(u{i},3)
%                rw{i}(j,k,l)=r{i}(j,k,l)*w{i}(j,k,l); 
%             end
%             
%             rwBar(j,k,i)=mean(rw{i}(j,k,:));
%         end
%     end
%     fprintf('Run %2.0f of rw done \n',i); 
% end
% 
% clear rw
% disp('Finished with rw'); 
% 
% for i =1:70
%     for j = 1:size(u{i},1)
%         for k = 1:size(u{i},2)
%             for l=1:size(u{i},3)
%                tw{i}(j,k,l)=t{i}(j,k,l)*w{i}(j,k,l); 
%             end
%             
%             twBar(j,k,i)=mean(tw{i}(j,k,:));
%         end
%     end
%     fprintf('Run %2.0f of tw done \n',i); 
% end
% 
% clear tw
% disp('Finished with tw');