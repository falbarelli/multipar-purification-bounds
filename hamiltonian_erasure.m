eta=.2;
bound_offdiag= @(d) (d-1).*(eta/(1-eta));
bound_diag= @(d) 4*((d-1)./d)*(eta/(1-eta));

dimvec=2:6; % rather slow for dim=6
bounds=[];
boundsReal=[];
boundsImag=[];
boundsDiag=[];

for dim=dimvec
    
    fprintf("Dimension %d \n", dim)
    
    dimIn=dim;
    dimOut=dim+1;

    npar=dim^2;
    id=eye(dim);
    
    %%% Kraus operators for the erasure channel
    nKraus=dim + 1;
    evec=[zeros(dim,1);1];
    KrausOps = zeros( dimOut * nKraus , dimIn );
    for i=1:dim
        KrausOps( (i-1)*dimOut +1 : i*dimOut , :) = sqrt(1-eta)*evec*id(i,:); 
    end
    % the last one is different
    KrausOps( dim*nKraus +1  : dimOut*nKraus, : )=sqrt(eta)*[id;zeros(1,dim)];
    
    KrausOpsDerivReal=zeros( dimOut*nKraus*(dim*(dim-1))/2 ,dimIn);
    p=0;
    for i=1:dim % real gens
        for j=i+1:dim
            p=p+1;
            for n=1:nKraus       
                KrausOpsDerivReal( (p-1)*nKraus*dimOut + (n-1)*dimOut + 1 : (p-1)*nKraus*dimOut + n*dimOut  ,  : ) =   -1j*KrausOps( (n-1)*dimOut + 1 : n*dimOut , :  )*GenGellMann(i-1,j-1,dim)/2;
            end
        end
    end
    
    KrausOpsDerivImag=zeros( dimOut*nKraus*(dim*(dim-1))/2 ,dimIn);
    p=0;
    for j=1:dim % imag gens
        for i=j+1:dim
            p=p+1;
            for n=1:nKraus       
                KrausOpsDerivImag( (p-1)*nKraus*dimOut + (n-1)*dimOut + 1 : (p-1)*nKraus*dimOut + n*dimOut  ,  : ) =   -1j*KrausOps( (n-1)*dimOut + 1 : n*dimOut , :  )*GenGellMann(i-1,j-1,dim)/2;
            end
        end
    end
    
    KrausOpsDerivDiag=zeros( dimOut*nKraus*dim ,dimIn);
    p=0;
    for j=1:dim % diag gens 
        p=p+1;
        for n=1:nKraus       
            KrausOpsDerivDiag( (p-1)*nKraus*dimOut + (n-1)*dimOut + 1 : (p-1)*nKraus*dimOut + n*dimOut  ,  : ) =   -1j*KrausOps( (n-1)*dimOut + 1 : n*dimOut , :  )*id(:,j)*id(j,:);
        end
    end
    
    boundsReal(end+1)=totalQFI_SDP(KrausOps,KrausOpsDerivReal,(dim*(dim-1))/2,nKraus,'asymptotic',dimIn,dimOut);
    boundsImag(end+1)=totalQFI_SDP(KrausOps,KrausOpsDerivImag,(dim*(dim-1))/2,nKraus,'asymptotic',dimIn,dimOut);
    boundsDiag(end+1)=totalQFI_SDP(KrausOps,KrausOpsDerivDiag,dim,nKraus,'asymptotic',dimIn,dimOut);
    
    KrausOpsDeriv = [KrausOpsDerivReal ; KrausOpsDerivImag ; KrausOpsDerivDiag ]; 
    bounds(end+1)=totalQFI_SDP(KrausOps,KrausOpsDeriv,npar,nKraus,'asymptotic',dimIn,dimOut);
    
end

%%% Confirmation of the results: these are all basically zero
bounds - (boundsReal+boundsImag+boundsDiag)
boundsReal-bound_offdiag(dimvec)
boundsImag-bound_offdiag(dimvec)
boundsDiag-bound_diag(dimvec)
