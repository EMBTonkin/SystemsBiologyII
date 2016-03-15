function [A, A1]=adjacency_pub(ncell,Perc_VIP,Perc_GABA, p_rand, V, G, cxn_scale, Self_cxn)
%%
A=zeros(ncell);
kk=1:ncell;
X=kk(ones(1,ncell),:);
s=sqrt(ncell);
kk=kk';
Var1=fix((kk-1)/s);
Var1=Var1(:,ones(1,ncell));
Var2=rem((kk-1),s);
Var2=Var2(:,ones(1,ncell));
A(1:ncell,1:ncell) = (sqrt((Var1-fix((X-1)/s)).^2+(Var2-rem((X-1),s)).^2)).\1;
A(A~=1)=0;
clear Var1 Var2 

c=rand(ncell); %%%creating symmetrical matrix of random numbers
for i=1:ncell
    for j=1:ncell
        c(i,j)=c(j,i);
    end
end

%%%%%%%%Generate random connectivity from Probability Density Fxn from Freeman (Neuron), 2013 %%%%%%%%%%%%%%%%%%
p=1:-log(1/(ncell/103*20.5))/0.176*ncell/103*cxn_scale; %%%the max number of connections allowed
q=20.5*exp(-0.175*103/ncell*p/cxn_scale); %%%First-order exponential distribution fit from Freeman data, scaled to ncell
pdf=q/sum(q);


for qq=2:length(pdf)  %%%makes probability density function
    pdf(qq)=pdf(qq)+pdf(qq-1);
end

cxn=ones(1,ncell);

for pp=1:ncell
    for ss=1:length(pdf)-1; 
        if p_rand(pp)>pdf(ss)
            cxn(pp)=ss+1;
        end
    end
end

bita_AVG=1-sqrt(1-mean(cxn(1:ncell))/ncell);
bita=1-(1-cxn(1:ncell)/ncell)/(1-bita_AVG);

%%%%%%% NJK: ADD CONNECTIONS RANDOMLY %%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% This method fails to generate enough cells with very low number of
%%% connections.  This is corrected in the loop after this one.

A=zeros(ncell);
for i=1:ncell
    for j=1:ncell
        if 1-(1-bita(i))*(1-bita(j))>c(i,j)
            A(i,j)=1;
        else A(i,j)=0;
        end
    end
end

%Fix distribution: assures that # of connections closely matches the pdf
%The number of iterations was sufficient at ncell=400
for iterations=1:25
    for adj=1:ncell
        while cxn(adj)<sum(A(:,adj))
            cut=fix(rand*ncell)+1;
            A(adj,cut)=0;
            A(cut,adj)=0;
        end
    end

    for adj=1:ncell
        while cxn(adj)>sum(A(:,adj))
            cut=fix(rand*ncell)+1;
            A(adj,cut)=1;
            A(cut,adj)=1;
        end
    end

%%%%%% EXCLUDE SELF CONNECTIONS (IF DESIRED)
    if Self_cxn==0    
        for i=1:ncell
            A(i,i)=0;
        end    
    end
end

%%%A(isnan(A))=0; % REPLACE NaN WITH 0

%%%%%%%%%%% Eliminating VIP connections based on which cells make VIP %%%%%
kk=1:ncell;
A1=A;
clear r k i w X 
% VIP PERCENTAGE
V(V<Perc_VIP)=0; %Perc_VIP = fraction (and not percentage) of VIP produced.
V(V~=0)=1;
Va=kk(V==1);
Vb=kk(V==0);
% GABA PERCENTAGE
G=G(1:ncell);
G(G<Perc_GABA)=0;
G(G~=0)=1;
Ga=kk(G==1);

for o = 1:length(Va)
    A(:,Va(o))=0;
end 
for o = 1:length(Ga)
    A1(:,Ga(o))=0;
end

%%% Making sure each cell always has at least one VIP connection.
%%% Giving connection to numerically closest VIP secreting neighbor
for i=1:ncell
    if sum(A(i,:),2)==0
        A(i,Vb(fix(rand*length(Vb))+1))=1;    
        while A(i,i)==1
            A(i,i)=0;
            A(i,Vb(fix(rand*length(Vb))+1))=1;
        end
    end
end    
clear Perc_VIP b wa s kk