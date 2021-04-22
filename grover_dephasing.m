eta=.7;
globBound = @(d ) ((4.*eta)./(1-eta)).*((d-1)./(d+2/eta));

dimvec=2:1:6;
ndim=length(dimvec);

bounds=zeros(ndim,1);
hmats=cell(ndim);

for k=1:ndim
    
    dim=dimvec(k);
    dimIn=dim;
    dimOut=dim;
    npar=dim;
    id=eye(dim);
    
    fprintf("Dimension %d \n", dim)

%%%%  Kraus operators used in the paper (non-minimal)
%     nKraus=dim+1;
%     KrausOps=zeros(dim*nKraus,dim);
%     KrausOps( 1 : dim , :) = sqrt( eta ) * eye(dim);
%     for i=2:nKraus
%         KrausOps( (i-1)*dim +1 : i*dim , :)=sqrt(1-eta)*id(:,i-1)*id(i-1,:);
%     end
    
%%%%  Alternatively we can use the canonical Kraus operators (obtained from
%%%%  diagonalizing the CJ matrix)
    DepCh=DephasingChannel(dim,eta);
    KrausDepCh=KrausOperators(DepCh);
    nKraus=dim;
    KrausOps=zeros(dim*nKraus,dim);
    for i=1:nKraus
        KrausOps( (i-1)*dim +1 : i*dim , :)=full(KrausDepCh{i});
    end
    
    KrausOpsDeriv = zeros( dimOut*nKraus*npar , dimIn );
    for p=1:npar
        for n=1:nKraus       
            KrausOpsDeriv( (p-1)*nKraus*dimOut + (n-1)*dimOut + 1 : (p-1)*nKraus*dimOut + n*dimOut  ,  : ) =   -1j*KrausOps( (n-1)*dimOut + 1 : n*dimOut , :  )*id(:,p)*id(p,:);
        end
    end
    
    [ bounds(k) , hmats{k}] = totalQFI_SDP(KrausOps,KrausOpsDeriv,npar,nKraus,'asymptotic',dimIn,dimOut);
end


plot(dimvec,bounds,'-o',dimvec,globBound(dimvec),'-x');