%%% This script reproduces the simulations used for Figure 3 and Supplementary 

%%%%%%%% Basics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Figures 4 and 5.
%%% Supplemental Figure 5A and B can also be reproduced in this script by
%%% modulating v_VIP through the parameter "rows".  NJK 1/21/15

%%% This script has paramter Shift_start adjusted by 12h and Phase_start is
%%% set to -9 rather than 12 to interrupt a 12h light phase by 9 hours and 
%%% test the dark pulse VRC prediction that re-entrainment times would be
%%% universally greater at this setting.  NJK 12/8/15

clear all; close all;
load hetparam.mat
xaa=5; %%% # of consecutive simulations; set to 1 to run once
col=[0,1];
rows=0.5;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     ;
iii=length(rows)*length(col); %%% # of runs per simulation;
End_of_simulation=1500; 
ncell = 400;
%%%%%%% Pre-allocating space for some outputs %%%%%%%%%%%%%%%%
scale=zeros(1,xaa*iii);
Mean_connections=zeros(1,xaa*iii);
Non_oscil=zeros(1,xaa*iii);
SI=zeros(1,xaa*iii);
Av=zeros(1,xaa*iii);
MaxMP=zeros(1,xaa*iii);
MinMP=zeros(1,xaa*iii);        
Maxv=zeros(1,xaa*iii);
Minv=zeros(1,xaa*iii);
STD=zeros(1,xaa*iii); 
SIvt=zeros(2*xaa*iii,round(End_of_simulation/24+1)); %%%For quickly looking at the sync index over time for each simulation

Sz=nan(ncell,iii*xaa);
ampz=nan(ncell,iii*xaa);
vz=nan(End_of_simulation*10+1,ncell,iii);
MPz=nan(End_of_simulation*10+1,ncell,iii);
%%
w=1;
while w<=xaa  
    clearvars -except rows col iii scale Mean_connections SI Av MaxMP Maxv Minv STD vsP0 vsB vmB vsP0z vsBz vmBz amp_intr p_intr Az A1z SIvt Sz ampz phdiff SI_mean av_amp av_S max_v vz MPz w End_of_simulation ncell xaa delay resync A A1
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);
    Perc_VIP = 0.2;
    Perc_GABA = 1; 
    cxn_scale=2.5; 
    Self_cxn=1;
    
    %Timing
    Coupling_start=150;
    Light_start=330; %%%When does the light swtich start (after coupling)
    light_hours_1=12; %%% hours of light per 24h day, before phase shift
    light_hours_2=12; %%% hours of light per 24h day, after phase shift
    Shift_start=690+12; %%% How long after light entrainment to perform phase shift
    Phase_shift=-9;  %%% How many hours to advance earlier the 24h circadian rhythm. Positive shifts lengthen the intermediate periods (phase delay).
    


    %%%SIMULATION
    ns=21;    
    p_rand=rand(1,ncell); %%%Assigns which cells will have the most connections
    V_rand=rand(1,ncell); %%%Assigns which cells will be VIP producers
    G_rand=rand(1,ncell); %%%Assigns which cells will be GABA producers
    l_rand=ones(ncell,1); %%%Assigns which cells will be sensitive to light
    sd = 0.01; 
%     vsP0= 0.94+sqrt((0.94*sd*6)^2)*randn(ncell,1);
%     vsB = 1.0*ones(ncell,1)+sqrt((1.0*sd*1)^2)*randn(ncell,1);
%     vmB = 0.8*ones(ncell,1)+sqrt((0.8*sd*1)^2)*randn(ncell,1);
    ba=0;
    vvip=0;
    xa=1;
    Cl_o=1;
    Cl_heterogeneity=0;
    Cl_distribution=Cl_o+Cl_heterogeneity*sd*randn(ncell,1);

    %%%Adjacency matrix for multicellular model
    [A,A1]=adjacency_pub(ncell,Perc_VIP,Perc_GABA,p_rand,V_rand,G_rand,cxn_scale,Self_cxn);% A--> VIP and A1-->GABA
    sumal=1./sum(A,2)'; % makes sure that there is no NaN in g(i) in goldbeter file
    sumal(isinf(sumal))=0;
    sumalGABA=1./sum(A1,2)';
    sumalGABA(isinf(sumalGABA))=0;
    fail=0;
    %%
    for ii=1:iii        
        v_VIP=rows(rem((ii-1),length(rows))+1);
        GABA_switch=col(fix((ii-1)/length(rows))+1);
        GABA_background=10^(-10)+GABA_switch*0.1;
        v_GABA=GABA_switch/2;
        t = 0:0.1:End_of_simulation;
        Lp_start=inf;Lp_end=inf;
        %%
        %%%%%% Solve ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        clear AA AA1 BC CC Ca_in Cl_in Fs GABA IN MB MP PC PCP PCC PCCP PCN PCNP Q Q1 S ST S1 SS S_GABA S_VIP vVIP X1 X2 Y1 aa ans...
        c c1 c2 c3 gGABA period period1 temp2 tim time v y Cl_distribution_transposed
        tic
        Simulation=(iii*(w-1)+ii)
        [v_VIP,v_GABA]
        if fail==0
            [t,y]=ode23(@ODEs_pub,t,IC16_pub(ncell,ns,GABA_switch),options,P16_pub,xa,ba,GABA_background,A,sumal,...
            A1,sumalGABA,ncell,ns,vsP0,vsB,vmB,Cl_distribution,vvip, light_hours_1,light_hours_2, Coupling_start,...
            Light_start, Shift_start, Phase_shift,End_of_simulation,l_rand,v_VIP,v_GABA,Lp_start,Lp_end);
        end
        toc
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        %%%%%%%% Recall all unique terms solved by the ODE %%%%%%%%%%%%    
        Ca_in = zeros(size(y,1),ncell);
        Ca_store = zeros(size(y,1),ncell);
        MP = zeros(size(y,1),ncell);
        MC = zeros(size(y,1),ncell);
        MB = zeros(size(y,1),ncell);
        PC = zeros(size(y,1),ncell);
        CC = zeros(size(y,1),ncell);
        PCP = zeros(size(y,1),ncell);
        CCP = zeros(size(y,1),ncell);
        PCC = zeros(size(y,1),ncell);
        PCN = zeros(size(y,1),ncell);
        PCCP = zeros(size(y,1),ncell);
        PCNP = zeros(size(y,1),ncell);
        BC = zeros(size(y,1),ncell);
        BCP = zeros(size(y,1),ncell);
        BN = zeros(size(y,1),ncell);
        BNP = zeros(size(y,1),ncell);
        IN = zeros(size(y,1),ncell);
        CB = zeros(size(y,1),ncell);
        vVIP = zeros(size(y,1),ncell);
        gGABA = zeros(size(y,1),ncell);
        
        for j=1:ncell
           Ca_in(:,j)= y(:,(j-1)*ns+1);
           Ca_store(:,j)= y(:,(j-1)*ns+2);
           MP(:,j)= y(:,(j-1)*ns+3);
           MC(:,j)=y(:,(j-1)*ns+4);
           MB(:,j)=y(:,(j-1)*ns+5);
           PC(:,j)=y(:,(j-1)*ns+6);
           CC(:,j)=y(:,(j-1)*ns+7);
           PCP(:,j)= y(:,(j-1)*ns+8);
           CCP(:,j)= y(:,(j-1)*ns+9);
           PCC(:,j)= y(:,(j-1)*ns+10);
           PCN(:,j)= y(:,(j-1)*ns+11);
           PCCP(:,j)= y(:,(j-1)*ns+12);
           PCNP(:,j)= y(:,(j-1)*ns+13);
           BC(:,j)= y(:,(j-1)*ns+14);
           BCP(:,j)= y(:,(j-1)*ns+15);
           BN(:,j)= y(:,(j-1)*ns+16);
           BNP(:,j)= y(:,(j-1)*ns+17);
           IN(:,j)= y(:,(j-1)*ns+18);
           CB(:,j)=y(:,(j-1)*ns+19);
           vVIP(:,j)=y(:,(j-1)*ns+20);
           gGABA(:,j)=y(:,(j-1)*ns+21);
        end
        
        if length(MP)==End_of_simulation*10+1
           MPz(:,:,ii)=MP;
           vsP0z(:,iii*(w-1)+ii)=vsP0;
           vsBz(:,iii*(w-1)+ii)=vsB;
           vmBz(:,iii*(w-1)+ii)=vmB;
           Az(:,:,iii*(w-1)+ii)=A;
           A1z(:,:,iii*(w-1)+ii)=A1;
         %%   
           SIW(:,:,ii)=SyncIndexWavos(t,MP);
           resync(iii*(w-1)+ii)=find(SIW(702+24:end,3,ii)>mean(SIW(702-96:690,3,ii))*0.95,1)+24
%%
           for j = 1:ncell;
              %        calculate period
               Y1(:,j)=filter(ones(1,240)/240,1,MPz(:,j,ii));
               Q(:,j)=MPz(241:length(t),j,ii)-Y1(241:length(t),j);%%% detrend data sets and subtract 24 h
               Q1(:,j)=filter(ones(1,30)/30,1,Q(:,j)); %%%3h running average of detrended set
               AA1(:,j)=filter(ones(1,241)/241,1,Q(:,j)); %%%24 h running average of detrended set
           end
           AA=AA1-Q1;
           q=zeros(length(AA)-240,ncell);
        

%%%% PERIOD AND AMP CALCULATIONS
            clear ST st period time tim S mntmp mxtmp amp q time tim c c3 S
            amp=nan(ncell,1);
            S=nan(ncell,1);  
            for j=1:ncell
                q(:,j)=AA(241:length(AA),j).*AA((241:length(AA))-1,j);%%% subtract 24h for 2nd time
                time(j).t = (find(q(:,j)<=0))+238;%%% you add 24 h back in order for the dimensions to match 'AA' dimensions.
                tim(j).t=real(time(j).t+(AA(time(j).t,j))./(AA(time(j).t,j)-AA(time(j).t+1,j)));%%% interpolate
                for tt=2:length(tim(j).t)  %%% Protected against rare error
                    if tim(j).t(tt)<tim(j).t(tt-1)
                        tim(j).t(tt)=time(j).t(tt)+tim(j).t(tt-1)-time(j).t(tt-1);
                    end
                end
                while tim(j).t(end)>length(AA)  %%%Protected against rare error
                    tim(j).t=tim(j).t(1:end-1);
                end

            %%%%%%% Coupled population (three cycles after coupling)%%%%%%%%%%
            if End_of_simulation>Coupling_start+72
            c(j).t=find(tim(j).t>=Coupling_start*10+720)-1 ;
            else c(j).t=find(tim(j).t>=720);
            end
             
            if  (length(c(j).t)>=3) && (c(j).t(1)>0)
                c3(j).t=abs(AA(round((tim(j).t(c(j).t-1)+tim(j).t(c(j).t-2))/2),j)-AA(round((tim(j).t(c(j).t-1)+tim(j).t(c(j).t))/2),j)); %find amplitude.
            else
                c3(j).t=0;
            end      
            if  mean(c3(j).t)>0.05
                period(j).t=(tim(j).t(c(j).t)-tim(j).t(c(j).t-2))/10;
                S(j,1)= mean(period(j).t(end-5:end));
                Sstd=std(period(j).t(end-5:end));
%                 if Sstd>0.2
%                     S(j,1)=nan;
%                 end
            else
                period(j).t=0;
                S(j,1)=nan;
            end     

            
                %%%Amp for each cell
                if MP(round(tim(j).t(end))+66,j)>MP(round(tim(j).t(end-1))+20,j)
                    mntmp=MP(round(tim(j).t([end-1,end-3,end-5,end-7,end-9]))+20,j);
                    mxtmp=MP(round(tim(j).t([end,end-2,end-4,end-6,end-8]))+66,j);
                else
                    mxtmp=MP(round(tim(j).t([end-1,end-3,end-5,end-7,end-9]))+66,j);
                    mntmp=MP(round(tim(j).t([end,end-2,end-4,end-6,end-8]))+20,j);
                end
                amp(j,1)=mean(mxtmp-mntmp);
                if std(mxtmp-mntmp)>0.05*amp(j)
                    amp(j)=nan;
                end        
            end     
            Average=mean(S(S>0));
            ss=std(S(S>0));
            bb=find(S~=0);
            aa=find(S==0);
            if length(MP)==End_of_simulation*10+1
                Sz(:,iii*(w-1)+ii)=S;
                ampz(:,iii*(w-1)+ii)=amp;
            end
%             for j=1:ncell
%                 timz(:,j,iii*(w-1)+ii).t=tim(j).t;
%                 periodz(:,j,iii*(w-1)+ii).t=period(j).t;
%             end
        
        %%%%%%%Synchronicity%%%%%%%%%%%%%%%%%%%%%% 
            temp2 = SyncIndex_pub(MP(100:length(t),bb)',t,mean(S(S>0)));
            for z=1:length(temp2)
                SIvt(2*iii*(w-1)+2*ii-1,z)=temp2(1,z); %%%SI vs. time: time component
                SIvt(2*iii*(w-1)+2*ii,z)=temp2(2,z); %%%SI vs. time: SI component
            end
        else fail=1;
        end
    end
    if fail==0
        %%
        %%%%%%%    ANALYSIS    %%%%%%%
     
        clear pks lks lksav jet p1 p2 ciptd mn mx phstart
        for ii=1:iii
            jet=jet;
            phstart=SIvt(2*(iii*(w-1)+ii)-1,SIvt(2*(iii*(w-1)+ii)-1,:)>Shift_start-72);
            ss=1:length(SIvt);
            lksav=(phstart(1)+(0:fix((End_of_simulation-SIvt(2*(iii*(w-1)+ii)-1,min(ss(SIvt(2*(iii*(w-1)+ii)-1,:)>Shift_start))-4))/24))*24)';
        end
        for ii=1:iii
            lks=nan(length(lksav),1);
            pks=nan(length(lksav),1);
            for j=1:ncell
                [pks,lks]=findpeaks(real(MPz(Shift_start*10-720:end,j,ii)));
                lks=lks(pks>1);
                lks=lks/10+Shift_start-72;
                dbl=find(diff(lks)<6);
                for d=1:length(dbl)
                    lks(dbl(d))=mean(lks(dbl(d))+lks(dbl(d)+1))/2;
                    lks(dbl(d)+1:end-1)=lks(dbl(d)+2:end);
                    lks=lks(1:end-1);
                    dbl=dbl-1;
                end
                phdiff(j,:,iii*(w-1)+ii)=lksav(1:end-3)-lks(1:length(lksav)-3);
                if phdiff(j,1,iii*(w-1)+ii)<-12
                    phdiff(j,:,iii*(w-1)+ii)=phdiff(j,:,iii*(w-1)+ii)+24;
                end
                if phdiff(j,1,iii*(w-1)+ii)>12
                    phdiff(j,:,iii*(w-1)+ii)=phdiff(j,:,iii*(w-1)+ii)-24;
                end
                for eee=1:length(phdiff(j,:,iii*(w-1)+ii))
                     while phdiff(j,eee,iii*(w-1)+ii)<-24
                         phdiff(j,eee,iii*(w-1)+ii)=phdiff(j,eee,iii*(w-1)+ii)+24;
                     end
                     while phdiff(j,eee,iii*(w-1)+ii)>24
                         phdiff(j,eee,iii*(w-1)+ii)=phdiff(j,eee,iii*(w-1)+ii)-24;
                     end
                end
           end
        end
        w=w+1;
    end
end
%%
if v_VIP<1
    index=v_VIP*10-4;
elseif v_VIP==1.5;
    index=6;
elseif v_VIP==3;
    index=7;
end
phdiff(rem(find(phdiff(:,4,1:iii)>4.2),400),:,:)=nan;
phdiff(rem(find(phdiff(:,4,1:iii)<-2.7),400),:,:)=nan;
kk=1:400;f=kk(isfinite(phdiff(:,end)))';
jack=nan(ncell,xaa*iii);
for ii=1:xaa*iii;
    for j=1:length(f)
        jack(f(j),ii)=find(abs(phdiff(f(j),5:end,ii))>abs(mean(phdiff(f(j),end-6:end,ii)))-0.5,1);
        if phdiff(f(j),end,ii)>6
            jack(f(j),ii)=jack(f(j),ii)-1;
        end
    end
    jackav(index,ii)=mean(jack(isfinite(jack(:,ii)),ii));
end
jackG(index,3:4)=[mean(jackav(2:2:iii)),std(jackav(2:2:iii))];
jackG(index,1:2)=[mean(jackav(1:2:iii-1)),std(jackav(1:2:iii-1))];
jackG(index,5)=ttest(jackav(2:2:iii),jackav(1:2:iii-1));

colr=['x','+'];
color=nan(ncell,1);
pdnoGav=mean(phdiff(:,:,1:2:xaa*length(col)-1),3);
pdwGav=mean(phdiff(:,:,2:2:xaa*length(col)),3);
pdnoGstd=std(phdiff(:,:,1:2:xaa*length(col)-1),[],3);
pdwGstd=std(phdiff(:,:,2:2:xaa*length(col)),[],3);
kk=1:400;
% pdnoGav(kk(pdnoGstd(:,end)>1),end)=nan;
% pdwGav(kk(pdwGstd(:,end)>1),end)=nan;
% f=intersect(kk(isfinite(pdnoGav(:,end))),kk(isfinite(pdwGav(:,end))));
color=nan(ncell,2);
for j=1:length(f)
    if pdnoGav(f(j),end)>0
        color(f(j),1)=1;
    elseif pdnoGav(f(j),end)<0
        color(f(j),1)=2;
    else color(f(j),1)=nan;
    end
    
    if pdwGav(f(j),end)>0
        color(f(j),2)=1;
    elseif pdwGav(f(j),end)<0
        color(f(j),2)=2;
    else color(f(j),2)=nan;
    end
end
jet=jet;
jac=nan(ncell,4);
for j=1:length(f)
    jac(f(j),1)=mean(jack(f(j),1:2:iii-1));
    jac(f(j),2)=std(jack(f(j),1:2:iii-1));
    jac(f(j),3)=mean(jack(f(j),2:2:iii));
    jac(f(j),4)=std(jack(f(j),2:2:iii));
end
mn=min(min([jac(f,1),jac(f,3)]));
mx=max(max([jac(f,1),jac(f,3)]));
kk=1:400;
p1=nan(ncell,6);

%%%%%%%%%%SUPPLEMENTARY FIGURE 4B%%%%%%%%%%%%%%
f=kk(isfinite(jac(:,1))==1);
ff=kk(isfinite(jac(:,1))==0);
for j=1:length(f)
    p1(f(j),1)=round((length(jet)-1)*((jac(f(j),1)-mn)/(mx-mn)))+1;
    figure(10000+v_VIP*100);hold on; plot(amp_intr(f(j)),p_intr(f(j)),'Color',jet(p1(f(j),1),:),'Marker',colr(color(f(j),1)),'LineStyle','none');
end
for j=1:length(ff)
    figure(10000+v_VIP*100);hold on; plot(amp_intr(ff(j)),p_intr(ff(j)),'Color','k','Marker','*','LineStyle','none');
end
set(gca,'YLim',[22,25],'XLim',[0,4.5])
xlabel('Intrinsic {\itPer} mRNA amplitude (nM)','FontSize',16)
ylabel('Intrinsic period (h)','FontSize',16)
%title(['v_V_I_P= ',num2str(v_VIP),' nM*h^-^1, v_G_A_B_A= 0 nM*h^-^1'],'FontSize',16)
colormap(jet);
cbarmax=colorbar('YTick',[1,length(jet)/2,length(jet)],'TickLabels',{round(mn*100)/100,round((mn+mx)/2*100)/100,round(mx*100)/100},'FontSize',16);
set(get(cbarmax,'ylabel'),'string','Time to Re-Entrain (d)','FontSize',16);
xlim=xlim;ylim=ylim;
line([median(amp_intr(isfinite(amp_intr))),median(amp_intr(isfinite(amp_intr)))],[ylim(1),ylim(2)],'LineStyle','--','Color','k');
line([xlim(1),xlim(2)],[median(p_intr(isfinite(p_intr))),median(p_intr(isfinite(p_intr)))],'LineStyle','--','Color','k');

%%%%%%%%%%SUPPLEMENTARY FIGURE 4C%%%%%%%%%%%%%%
f=kk(isfinite(jac(:,3))==1);
ff=kk(isfinite(jac(:,3))==0);
for j=1:length(f)
    p1(f(j),2)=round((length(jet)-1)*((jac(f(j),3)-mn)/(mx-mn)))+1;
    figure(20000+v_VIP*100);hold on; plot(amp_intr(f(j)),p_intr(f(j)),'Color',jet(p1(f(j),2),:),'Marker',colr(color(f(j),2)),'LineStyle','none');
end
for j=1:length(ff)
    figure(20000+v_VIP*100);hold on; plot(amp_intr(ff(j)),p_intr(ff(j)),'Color','k','Marker','*','LineStyle','none');
end
set(gca,'YLim',[22,25],'XLim',[0,4.5])
xlabel('Intrinsic {\itPer} mRNA amplitude (nM)','FontSize',16)
ylabel('Intrinsic period (h)','FontSize',16)
%title(['v_V_I_P= ',num2str(v_VIP),' nM*h^-^1, v_G_A_B_A= 0.5 nM*h^-^1'],'FontSize',16)
colormap(jet);
cbarmax=colorbar('YTick',[1,length(jet)/2,length(jet)],'TickLabels',{round(mn*100)/100,round((mn+mx)/2*100)/100,round(mx*100)/100},'FontSize',16);
set(get(cbarmax,'ylabel'),'string','Time to Re-Entrain (d)','FontSize',16);
line([median(amp_intr(isfinite(amp_intr))),median(amp_intr(isfinite(amp_intr)))],[ylim(1),ylim(2)],'LineStyle','--','Color','k');
line([xlim(1),xlim(2)],[median(p_intr(isfinite(p_intr))),median(p_intr(isfinite(p_intr)))],'LineStyle','--','Color','k');

%%%%%%%%%%SUPPLEMENTARY FIGURE 4D%%%%%%%%%%%%%%
f=kk(isfinite(jac(:,1)));
ff=kk(isfinite(jac(:,1))==0);
mn=-max([abs(min(-jac(f,1)+jac(f,3))),abs(max(-jac(f,1)+jac(f,3)))]);
mx=max([abs(min(-jac(f,1)+jac(f,3))),abs(max(-jac(f,1)+jac(f,3)))]);
jdiff=nan(ncell,1);
for j=1:length(f)
    jdiff(f(j))=-jac(f(j),1)+jac(f(j),3);
    if jdiff(f(j))>=0
        p1(f(j),3)=round((length(jet)-1)*(jdiff(f(j))/mx/2+0.5))+1;
    elseif jdiff(f(j))<0
        p1(f(j),3)=round(-(jdiff(f(j))-mn)/mn/2*(length(jet)-1))+1;
    end
end
for j=1:length(f)
    figure(30000+v_VIP*100);hold on; plot(amp_intr(f(j)),p_intr(f(j)),'Color',jet(p1(f(j),3),:),'Marker','.','LineStyle','none');
end
for j=1:length(ff)
    figure(30000+v_VIP*100);hold on; plot(amp_intr(ff(j)),p_intr(ff(j)),'Color','k','Marker','*','LineStyle','none');
end
set(gca,'YLim',[22,25],'XLim',[0,4.5])
xlabel('Intrinsic {\itPer} mRNA amplitude (nM)','FontSize',16)
ylabel('Intrinsic period (h)','FontSize',16)
%title(['Difference with GABA: v_V_I_P= ',num2str(v_VIP),' nM*h^-^1'],'FontSize',16)
colormap(jet);
cbarmax=colorbar('YTick',[1,length(jet)/2,length(jet)],'TickLabels',{round(mn*100)/100,0,round(mx*100)/100},'FontSize',16);
set(get(cbarmax,'ylabel'),'string','Difference in Time to Re-Entrain (d)','FontSize',16); 
line([median(amp_intr(isfinite(amp_intr))),median(amp_intr(isfinite(amp_intr)))],[ylim(1),ylim(2)],'LineStyle','--','Color','k');
line([xlim(1),xlim(2)],[median(p_intr(isfinite(p_intr))),median(p_intr(isfinite(p_intr)))],'LineStyle','--','Color','k');

%%%%%%%%%%SUPPLEMENTARY FIGURE 4A%%%%%%%%%%%%%%
colr=['r','b','g','y','w'];
for j=1:length(f)
    if pdnoGav(f(j),end)>0
        if pdwGav(f(j),end)>0
            color(f(j))=1;
        elseif pdwGav(f(j),end)<0
            color(f(j))=3;
        else color(f(j))=5;
        end
    elseif pdnoGav(f(j),end)<0
        if pdwGav(f(j),end)<0
            color(f(j))=2;
        elseif pdwGav(f(j),end)>0
            color(f(j))=4;
        else color(f(j))=5;
        end
    end
end
f=kk(isfinite(color(:,1))==1);
ff=kk(isfinite(color(:,1))==0);
for j=1:length(f)
    figure(100+v_VIP*100);hold on; plot(amp_intr(f(j)),p_intr(f(j)),'Color',colr(color(f(j))),'Marker','.','LineStyle','none');
end
for j=1:length(ff)
    figure(100+v_VIP*100);hold on; plot(amp_intr(ff(j)),p_intr(ff(j)),'Color','k','Marker','*','LineStyle','none');
end
set(gca,'YLim',[22,25],'XLim',[0,4.5])
xlabel('Intrinsic {\itPer} mRNA amplitude (nM)','FontSize',16)
ylabel('Intrinsic period (h)','FontSize',16)
%title(['v_V_I_P= ',num2str(v_VIP),'nM*h^-^1'],'FontSize',16)
line([median(amp_intr(isfinite(amp_intr))),median(amp_intr(isfinite(amp_intr)))],[ylim(1),ylim(2)],'LineStyle','--','Color','k');
line([xlim(1),xlim(2)],[median(p_intr(isfinite(p_intr))),median(p_intr(isfinite(p_intr)))],'LineStyle','--','Color','k');
% %%
% %%%%%%Additional tabulation for reporting numeric data within the text is
% %%%%%%performed here.
% jad(index,1:2)=[mean(mean(jack(find(color==1),1:2:iii-1))),std(mean(jack(find(color==1),1:2:iii-1)))];
% jad(index,3:4)=[mean(mean(jack(find(color==1),2:2:iii))),std(mean(jack(find(color==1),2:2:iii)))];
% jad(index,5:6)=[length(find(color==1)),ttest(mean(jack(find(color==1),2:2:iii)),mean(jack(find(color==1),1:2:iii-1)))];
% jad(index,7:8)=[mean(mean(jack(find(color==2),1:2:iii-1))),std(mean(jack(find(color==2),1:2:iii-1)))];
% jad(index,9:10)=[mean(mean(jack(find(color==2),2:2:iii))),std(mean(jack(find(color==2),2:2:iii)))];
% jad(index,11:12)=[length(find(color==2)),ttest(mean(jack(color==2,2:2:iii)),mean(jack(find(color==2),1:2:iii-1)))];
% jad(index,13:14)=[mean(mean(jack(find(color==3),1:2:iii-1))),std(mean(jack(find(color==3),1:2:iii-1)))];
% jad(index,15:16)=[mean(mean(jack(find(color==3),2:2:iii))),std(mean(jack(find(color==3),2:2:iii)))];
% jad(index,17:18)=[length(find(color==3)),ttest(mean(jack(find(color==3),2:2:iii)),mean(jack(find(color==3),1:2:iii-1)))];
% jad(index,19:20)=[mean(mean(jack(f,1:2:iii-1))),std(mean(jack(f,1:2:iii-1)))];
% jad(index,21:22)=[mean(mean(jack(f,2:2:iii))),std(mean(jack(f,2:2:iii)))];
% jad(index,23:24)=[length(f),ttest(mean(jack(f,1:2:iii-1)),mean(jack(f,2:2:iii)))];
% 
% for w=1:xaa
%     jackdiffncell(:,w)=jack(:,(w-1)*iii+1)-jack(:,w*iii);
% end
% jackdiff(index,1:4)=[mean(mean(jackdiffncell(find(color==1),:))),std(mean(jackdiffncell(find(color==1),:))),length(find(color==1)),ttest(mean(jackdiffncell(find(color==1),:),2))];
% jackdiff(index,9:12)=[mean(mean(jackdiffncell(find(color==3),:))),std(mean(jackdiffncell(find(color==3),:))),length(find(color==3)),ttest(mean(jackdiffncell(find(color==3),:),2))];
% jackdiff(index,5:8)=[mean(mean(jackdiffncell(find(color==2),:))),std(mean(jackdiffncell(find(color==2),:))),length(find(color==2)),ttest(mean(jackdiffncell(find(color==2),:),2))];
% jackdiff(index,13:16)=[mean(mean(jackdiffncell(f,:))),std(mean(jackdiffncell(f,:))),length(f),ttest(mean(jackdiffncell(f,:),2))];
% %%
% index=1;
% peaktimenoG(:,index)=-pdnoGav(:,4)+mean(SIvt(1:4:17,29))-702;
% peaktimewG(:,index)=-pdwGav(:,4)+mean(SIvt(3:4:19,29))-702;
% %%
% ad=nan(400,iii);
% for w=1:iii
%     for j=1:ncell
%         if phdiff(j,end,w)>0
%             ad(j,w)=1;
%         elseif phdiff(j,end,w)<0
%             ad(j,w)=2;
%         end
%     end
% end
%%
% %%%%%Plotting daylight for pre-simulation testing daylight timing parameters
% Light_period=24;
% Shift_start=690+12;
% Phase_shift=-9;  %%% How many hours to advance earlier the 24h circadian rhythm. Positive shifts lengthen the intermediate periods (phase delay).
% for t=0.1:0.1:End_of_simulation
%     if  (t>=Light_start) && (t<=Shift_start+Phase_shift)
%         if (t-Light_start)/Light_period-fix((t-Light_start)/Light_period)-light_hours_1/Light_period>=0
%             daylight(round(t*10))=0;
%         else
%             daylight(round(t*10))=1;
%         end
%     elseif (t>Shift_start+Phase_shift)
%          if (t-Light_start-Phase_shift)/Light_period-fix((t-Light_start-Phase_shift)/Light_period)-light_hours_2/Light_period>=0
%              daylight(round(t*10))=0;
%          else
%              daylight(round(t*10))=1;
%          end
%     else
%     daylight(round(t*10))=0;
%     end
% end
% figure(2);area(daylight)