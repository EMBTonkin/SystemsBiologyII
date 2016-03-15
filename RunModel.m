 xaa=0
 for w=1:length(xaa)
for ii=1:1
tic
%%%%%%%% Basics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('RelTol',1e-3,'AbsTol',1e-6);
ncell = 400;%# of cells
ns = 21;%# of states
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load A
% load A1
% load vmB
% load vsB
% load vsP0
%%%%%%%% Create Heterogeneity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% Kingsbury, Taylor, Henson, pg 2, 1.2. Network Heterogeneity
 sd = 0.01;
vsP0 =0.94+sqrt((0.94*sd*6)^2)*randn(ncell,1);
vsB = (1.0*ones(ncell,1))+ sqrt((1.0*sd*1)^2)*randn(ncell,1);
vmB = (0.8*ones(ncell,1))+ sqrt((0.8*sd*1)^2)*randn(ncell,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%%%%  Connectivity Matrix %%%%%%%%%%%%%%%%%%%%%%%%%%
 Perc_VIP =0.2 ;%* %Perc_VIP is %VIP producersb=rand(1,ncell);
 bita=0.05
 [A A1]=adjacency(ncell,bita,Perc_VIP);% A--> VIP and A1-->GABA
sumal=1./sum(A,2)'; % I did this to make sure that I don't get NaN in g(i) in goldbeter file
sumal(isinf(sumal))=0;
scale = mean(sum(A,2))
sumalGABA=1./sum(A1,2)';
sumalGABA(isinf(sumalGABA))=0;
mean(sum(A1,2))
%%%^%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Solve ODEs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%exit=xaa(w)
xa=1;% Controls Calcium content
ba=0;% Controls GABA content
Cl_o=1;
vvip=0;
t1=xaa(w)
t = [0:0.1:450];
[t,y]=ode23(@ODEs,t,IC16(ncell,ns),options,P16,xa,ba,A,sumal,A1,sumalGABA,ncell,ns,vsP0,vsB,vmB,Cl_o,vvip,t1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc
for j =1:ncell
 MP(:,j)= y(:,(j-1)*ns+3);
end
% mean_MP=mean(MP,2);
%  t_max(w)=6001-240+find(mean_MP(6001-240:6001)==max(mean_MP(6001-240:6001)))+1
% %  tmax=5913;
% %  Dt(w)=tmax-t_max(w)
% figure(1);plot(t,mean(MP,2),'k')
% hold on
% end
%  end
 
% k = 1.4*10^(-23) ;
% T = 37 + 273.15; % Kelvin
% q = 1.6*10^(-19) ;% C
% Ca_out = 5;%;%mM\
% Kex2=1;
% nex2=-1;
% Vex2=3.5;%
% % tic 
for j =1:ncell
MP(:,j)= y(:,(j-1)*ns+3);
Ca_in(:,j)= y(:,(j-1)*ns+1);
CC(:,j)=y(:,(j-1)*ns+7);
BC(:,j)=y(:,(j-1)*ns+14);
vVIP(:,j)=y(:,(j-1)*ns+20);
gGABA(:,j)=y(:,(j-1)*ns+21);
MB(:,j)=y(:,(j-1)*ns+5);
end
% % toc
% % pp=P16;
% % KD = pp(55);
% % clear pp j
% % 
% % 
% % S_GABA=zeros(length(t),ncell);
% % 
for j=1:length(t)
% %     VIP=vVIP(j,:);
% %     VIP=VIP(ones(1,ncell),:);
% %     S_VIP(j,:)=sum(VIP.*A,2)'.*sumal+vvip;
% %     clear VIP
% %  
    GABA=0.1+gGABA(j,:);
    GABA=GABA(ones(1,ncell),:);
    S_GABA(j,:)=ba+sum(GABA.*A1,2).*sumalGABA'; % GABA binding to cell surface
    clear GABA
end
% % v=zeros(length(t),ncell);
% % 
for j=1:length(t)
[v(j,:)]= FiringRates(Ca_in(j,:),(S_GABA(j,:)),Fir,CC(j,:),BC(j,:),MP(j,:),Cl_o);
end
% % clear j i  CC BC vVIP MB GABA  ECa VEX_ca Perc_VIP 
% % 
% % 
% % 
% % %%%%%%%%%%%%%% Period calc.%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%% For more info on cross-over analysis%%%
% % %%%%%%%%%%% go to Abe et al 2002 %%%%%%%%%%%%%%%%%
  for j = 1:ncell;
%        calculate period
        Y1(:,j)=filter(ones(1,240)/240,1,MP(:,j));
        Q(:,j)=MP(241:length(t),j)-Y1(241:length(t),j); % detrend data sets and subtract 24 h
        Q1(:,j)=filter(ones(1,30)/30,1,Q(:,j)); % 3h running average of detrended set
        AA1(:,j)=filter(ones(1,241)/241,1,Q(:,j)); %24 h running average of detrended set
  end
AA=AA1-Q1;
size(AA)
q=zeros(length(AA)-240,ncell);
  for j=1:ncell     
        q(:,j)=AA(241:length(AA),j).*AA((241:length(AA))-1,j);% substract 24h for 2nd time
        time(j).t = (find(q(:,j)<=0))+238;% you add 24 h back in order for the dimensions to match 'AA' dimensions.
        tim(j).t=time(j).t+(AA(time(j).t,j))./(AA(time(j).t,j)-AA(time(j).t+1,j));% interpolate

        %%%%%%%%%%%%uncoupled population%%%%%%%%%%%%%%%%%%%%%%%%%%
       c1(j).t=find(tim(j).t >=800 & tim(j).t<=1260);
       if (length(c1(j).t)>=3) && (c1(j).t(1)>0) %&& (mean(MP(round(tim(j).t(c1(j).t)),j))>=0.05)
           c2(j).t=abs(MP(round((tim(j).t(c1(j).t-1)+tim(j).t(c1(j).t-2))/2),j)-MP(round((tim(j).t(c1(j).t-1)+tim(j).t(c1(j).t))/2),j)); %find amplitude.
       else
           c2(j).t=0;
       end
       if c2(j).t(end)==0 || (max(MP(1200:1500,j))-min(MP(1200:1500,j)))<=0.15%abs(c2(j).t(end)-c2(j).t(end-2)) <=0.05
           period1(j).t=0;
       else
           period1(j).t=(tim(j).t(c1(j).t)-tim(j).t(c1(j).t-2))/10;
       end
       S1(j,1)=mean(period1(j).t);
     
       Average1=mean(S1(S1>0));
       ss1=std(S1(S1>0));
       X1=find(S1==0);
       X2=find(S1>0);     
     
       %%%%%%% Coupled population (three cycles after coupling)%%%%%%%%%%
        c(j).t=find(tim(j).t>=2500)-1 ;
        if  (length(c(j).t)>=5) && (c(j).t(1)>0) %&&(mean(MP(round(tim(j).t(c(j).t)),j))>= 0.05)
            c3(j).t=abs(AA(round((tim(j).t(c(j).t-1)+tim(j).t(c(j).t-2))/2),j)-AA(round((tim(j).t(c(j).t-1)+tim(j).t(c(j).t))/2),j)); %find amplitude.
        else
            c3(j).t=0;
        end

        if  mean(c3(j).t)>0.05
            period(j).t=(tim(j).t(c(j).t)-tim(j).t(c(j).t-2))/10;
        else
            period(j).t=0;
        end
        S(j,1)= (period(j).t(end));
        Average=mean(S(S>0));
        ss=std(S(S>0));
  end
     bb=find(S~=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);plot(t,MP(:,bb))
aa=find(S==0);
% figure(2);plot(t,MP(:,aa))
Non_oscil(w,ii)=length(find(S==0))/ncell*100
  clear Y1 Q Q1 AA1  AA q time.t tim.t c3 c c1 c2 tim time  period1 period  ss1 

%%%%%%%Synchronicity%%%%%%%%%%%%%%%%%%%%%%
temp2 = SyncIndex(MP(100:length(t),bb)',t,mean(S(S>0)))
SI(w,ii)  = max(temp2(length(temp2)),temp2(length(temp2)-1))
Av(w,ii)=Average

MaxMP(w,ii)=mean(max(real(MP(3000:end,:))))
MinMP(w,ii)=mean(min(real(MP(3000:end,:))))

Maxv(w,ii)=mean(max(real(v(3000:end,:))))
Minv(w,ii)=mean(min(real(v(3000:end,:))))


SS(w,ii)=size(X1,1)/ncell
STD(w,ii)=std(S(S>0))

%%%%%%%%% Amplitude%%%%%%%%%%
% %     for j=1:length(bb)
% %        T0(j)=round(t(end)-S(bb(j)));
% %        minimumMP(j)=find(MP(T0(j)*10:end,bb(j))==min(MP(T0(j)*10:end,bb(j))));
% %        maximumMP(j)=find(MP(T0(j)*10:end,bb(j))==max(MP(T0(j)*10:end,bb(j))));
% %        amplitudeMP(j)= MP(T0(j)*10 + maximumMP(j)-1,bb(j)) - MP(T0(j)*10 + minimumMP(j)-1,bb(j));
% %        MAXMP(j)= MP(T0(j)*10 + maximumMP(j)-1,bb(j));
% %        MINMP(j)=MP(T0(j)*10 + minimumMP(j)-1,bb(j));    
% %        
% %        minimumv(j)=find(v(T0(j)*10:end,bb(j))==min(v(T0(j)*10:end,bb(j))));
% %        maximumv(j)=find(v(T0(j)*10:end,bb(j))==max(v(T0(j)*10:end,bb(j))));
% %        MAXV(j)=v(T0(j)*10 + maximumv(j)-1,bb(j));
% %        MINV(j)=v(T0(j)*10 + minimumv(j)-1,bb(j));
% %        amplitudev(j)= v(T0(j)*10 + maximumv(j)-1,bb(j)) - v(T0(j)*10 + minimumv(j)-1,bb(j));
% %        
% %        minimumvrest(j)=find(Vrest(T0(j)*10:end,bb(j))==min(Vrest(T0(j)*10:end,bb(j))));
% %        maximumvrest(j)=find(Vrest(T0(j)*10:end,bb(j))==max(Vrest(T0(j)*10:end,bb(j))));
% %        MAXVrest(j)=Vrest(T0(j)*10 + maximumvrest(j)-1,bb(j));
% %        MINVrest(j)=Vrest(T0(j)*10 + minimumvrest(j)-1,bb(j));
% %        amplitudevrest(j)= Vrest(T0(j)*10 + maximumvrest(j)-1,bb(j)) - Vrest(T0(j)*10 + minimumvrest(j)-1,bb(j));
% %     
% %        minimumvIP(j)=find(IPSP(T0(j)*10:end,bb(j))==min(IPSP(T0(j)*10:end,bb(j))));
% %        maximumvIP(j)=find(IPSP(T0(j)*10:end,bb(j))==max(IPSP(T0(j)*10:end,bb(j))));
% %        MAXIPSC(j)=IPSP(T0(j)*10 + maximumvIP(j)-1,bb(j));
% %        MINIPSC(j)=IPSP(T0(j)*10 + minimumvIP(j)-1,bb(j));
% %        amplitudeIPSC(j)= IPSP(T0(j)*10 + maximumvIP(j)-1,bb(j)) - IPSP(T0(j)*10 + minimumvIP(j)-1,bb(j));
% % 
% %     end
% %     
% %     Mean_minMP(w,ii)=mean(MINMP)
% %     Mean_maxMP(w,ii)=mean(MAXMP)
% %     Mean_minfiri(w,ii)=mean(MINV)
% %     Mean_maxfiri(w,ii)=mean(MAXV)
% %     Mean_minverst(w,ii)=mean(MINVrest)
% %     Mean_maxverst(w,ii)=mean(MAXVrest)    
% %     Mean_minIPSC(w,ii)=mean(MINIPSC)
% %     Mean_maxIPSC(w,ii)=mean(MAXIPSC)    
% % 
% %     clear minimumMP minimumv maximumMP maximumv T0 j
% %     
% % Mean_amplMP(w,ii)=mean(amplitudeMP)
% % Mean_amplv(w,ii)=mean(amplitudev)
% % Mean_amplvrest(w,ii)=mean(amplitudevrest)
% % Mean_amplIP(w,ii)=mean(amplitudeIPSC)
% % Mean_levels_Vrest(w,ii)=mean(mean(Vrest(3000:end,:)))
end
 end
% %  
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %   
% % figure(w);plot(t,vVIP)
% % 
% % 
% % 
% % clear Average1 ss1 X1 X2 y c c2 period 
% % clear MP Ca_in CC BC vVIP MB GABA S
% % figure(w);plot(t,vVIP)
% % 
% % figure(1);plot(t,MP)
% %  
% %  
% % % figure(1);plot(t,MP)
% %   figure(1);plot(t,mean(MP,2),'r')
% %   hold on
% %  figure(2);plot(t,mean(v,2),'r')
% %  hold on
% % 
% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
% % return
% % a1(1,1)=min(v(2300:2300+600));
% % a1(2,1)=max(v(2300:2300+600));
% % a1(3,1)=mean(v(2300:2300+600));
% % 
% % a1(1,2)=min(Ca_in(2300:2300+600));
% % a1(2,2)=max(Ca_in(2300:2300+600));
% % a1(3,2)=mean(Ca_in(2300:2300+600));
% % 
% % a1(1,3)=min(MP(2300:2300+600));
% % a1(2,3)=max(MP(2300:2300+600));
% % a1(3,3)=mean(MP(2300:2300+600));