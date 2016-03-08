function SI = SyncIndex(M,t,bin)
% OC = original comment
% equations refer to Vasalou 2011 A Multicellular Model For Differential regulation, page 5, "Computational studies and data analysis" section


%OC: ncell=123;
temp = size(M);
nbin = round(t(temp(2))/bin)-1;
ncell = temp(1);
TOL = 0.3;


aM=mean(M);

%OC: compute the phase angles
ind(1)=1;
for i = 1:nbin;        %OC: determines the reference
    temp = find(abs(t-i*bin)<TOL);
    ind(i+1) = temp(1);
    ind2 = find((aM(ind(i):ind(i+1)-1)==max(aM(ind(i):ind(i+1)-1))));
    tim(i) = t(ind2+ind(i));
    for j = 1:ncell;       
        ind3 = find(M(j,ind(i):ind(i+1)-1)==max(M(j,ind(i):ind(i+1)-1)));
        time(i,j) = t(ind3+ind(i));
        theta(i,j) = (time(i,j)-tim(i))/bin*2*pi;
    end;
end


%OC: compute the synchronization indices
for i = 1:nbin;
    OP(i) = exp(1i*theta(i,1));
    for j = 2:ncell;
        OP(i) = OP(i)+exp(1i*theta(i,j));
    end
    OP(i) = OP(i)/ncell;
% OC:   r(i) = abs(OP(i));
    day(i) = i;
end

SI = abs(OP);
