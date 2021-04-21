eta=.4;
BoundLoss = @(d) ((4*eta)/(1-eta))*(d - 1)./d;
IncompLoss = @(d) (d - 1).^2./(4*(d-2));

dimvec=2:1:6;
ndim=length(dimvec);
bounds=zeros(ndim,1);
incompat=zeros(ndim,1);
hmats=cell(ndim);

for k=1:ndim
    
    dim=dimvec(k);
    
    fprintf("Dimension %d \n", dim)
    
    dimIn=dim;
    dimOut=dim+1;
    
    npar=dim-1;
   
    id=eye(dim); 

    nKraus=dim+1;
    KrausOps=zeros(dimOut*nKraus,dimIn);
    KrausOps( 1 : dimOut , :) = sqrt( eta ) * [ eye(dim) ; zeros(1,dim) ];
    for i=2:nKraus
        KrausOps( (i-1)*dimOut +1 : i*dimOut , :)=sqrt(1-eta)*[ zeros(dim) ; id(i-1,:) ];
    end    
    
    wInv=zeros(npar,1);
    KrausOpsDeriv = zeros( dimOut*nKraus*npar , dimIn );
    for p=1:npar
        for n=1:nKraus       
            KrausOpsDeriv( (p-1)*nKraus*dimOut + (n-1)*dimOut + 1 : (p-1)*nKraus*dimOut + n*dimOut  ,  : ) =   -1j*KrausOps( (n-1)*dimOut + 1 : n*dimOut , :  )*id(:,p)*id(p,:);
        end
        wInv(p) = totalQFI_SDP(KrausOps,KrausOpsDeriv( (p-1)*nKraus*dimOut + 1 : p*nKraus*dimOut ,  : ),1,nKraus,'asymptotic',dimIn,dimOut);
    end
    
    [ bounds(k), hmats{k} ]   = totalQFI_SDP(KrausOps,KrausOpsDeriv,npar,nKraus,'asymptotic',dimIn,dimOut);
    incompat(k) = npar*(1/totalQFI_SDP(KrausOps,KrausOpsDeriv,npar,nKraus,'asymptotic',dimIn,dimOut,wInv));
end

%%% Confirmation of the results: these are all basically zero
bounds-BoundLoss(dimvec-1)  % the case dim=2 should be handled separately
incompat-IncompLoss(dimvec)
