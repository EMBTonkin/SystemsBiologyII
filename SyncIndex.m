function SI = SyncIndex_pub(M,t,bin)
clear timmit
temp = size(M);
nbin = round(t(temp(2))/bin)-1;
ncell = temp(1);
TOL = 0.3;


aM=mean(M);

%compute the phase angles
ind(1)=1;
for i = 1:nbin;        %determines the reference
    temp = find(abs(t-i*bin)<TOL);
    ind(i+1) = temp(1);
    ind2 = find((aM(ind(i):ind(i+1)-1)==max(aM(ind(i):ind(i+1)-1))));
    timmit(i) = t(ind2+ind(i));
    for j = 1:ncell;       
        ind3 = find(M(j,ind(i):ind(i+1)-1)==max(M(j,ind(i):ind(i+1)-1)));
        timemit(i,j) = t(ind3+ind(i));
        theta(i,j) = (timemit(i,j)-timmit(i))/bin*2*pi;
    end;
end

%compute the synchronization indices
for i = 1:nbin;
    OP(i) = exp(1i*theta(i,1));
    for j = 1:ncell;
        OP(i) = OP(i)+exp(1i*theta(i,j));
    end
    OP(i) = OP(i)/ncell;
%    r(i) = abs(OP(i));
    day(i) = i;
end

%Output shows the time of peak Per mRNA of the mean population as the first
%row and the corresponding syncronization index as the second row
SI(1,:) = timmit+10;
SI(2,:) = abs(OP);
