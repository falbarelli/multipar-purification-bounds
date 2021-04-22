sx=full(Pauli(1));
sy=full(Pauli(2));
sz=full(Pauli(3));

npar=2;
nKraus=4;
dimIn=2;
dimOut=dimIn;

NN=.25;

gammavals=0.05:.05:.95;
ngamma=length(gammavals);

bounds=zeros(ngamma,1);
boundsSing1=zeros(ngamma,1);
boundsSing2=zeros(ngamma,1);
pIncCost=zeros(ngamma,1);

hmats=cell(ngamma);
hmats2=cell(ngamma);
hmat1=cell(ngamma);
hmat2=cell(ngamma);

% RLD bound (valid for NN <= 1/2)
RLDbound = @(gamma) (NN+(-1).*NN.^2).^(-1)+(1/4).*gamma.^(-2).*((-4)+(1+(-1).*NN).^(-1)+(NN+(-1).*gamma.*NN).^(-1));

for k=1:ngamma
    
    gamma=gammavals(k);
    
    K0=sqrt(1-NN)*[1 , 0 ; 0 , sqrt(1-gamma)];
    K1=sqrt(1-NN)*[0, sqrt(gamma) ; 0 , 0];
    K2=sqrt(NN)*[sqrt(1-gamma),0;0,1];
    K3=sqrt(NN)*[0,0;sqrt(gamma) ,0 ];
    KrausOps=[ K0 ; K1 ; K2 ; K3 ];

    d1K0=-(0.5/sqrt(1-NN))*[1 , 0 ; 0 , sqrt(1-gamma)];
    d1K1=-(0.5/sqrt(1-NN))*[0, sqrt(gamma) ; 0 , 0];
    d1K2=(0.5/sqrt(NN))*[sqrt(1-gamma),0;0,1];
    d1K3=(0.5/sqrt(NN))*[0,0;sqrt(gamma) ,0 ];
    d1KrausOps=[ d1K0 ; d1K1 ; d1K2 ; d1K3 ];

    d2K0=sqrt(1-NN)*[0 , 0 ; 0 , -0.5/sqrt(1-gamma)];
    d2K1=sqrt(1-NN)*[0, 0.5/sqrt(gamma) ; 0 , 0];
    d2K2=sqrt(NN)*[-0.5/sqrt(1-gamma),0;0,0];
    d2K3=sqrt(NN)*[0,0; 0.5/sqrt(gamma) ,0 ];
    d2KrausOps=[ d2K0 ; d2K1 ; d2K2 ; d2K3 ];

    dNewKrausOps = [d1KrausOps;d2KrausOps];
    
    [ boundsSing1(k), hmat1{k}] = totalQFI_SDP(KrausOps,d1KrausOps,1,nKraus,'single',dimIn,dimOut);
    [ boundsSing2(k), hmat2{k}] = totalQFI_SDP(KrausOps,d2KrausOps,1,nKraus,'single',dimIn,dimOut);
        
    [ bounds(k), hmats{k}] = totalQFI_SDP(KrausOps,dNewKrausOps,npar,nKraus,'single',dimIn,dimOut);
    [pIncCost(k), hmats2{k}]= totalQFI_SDP(KrausOps,dNewKrausOps,npar,nKraus,'single',dimIn,dimOut , [boundsSing1(k),boundsSing2(k)] );
    pIncCost(k)=2/pIncCost(k);

end

figure
plot(gammavals,bounds,'-x',gammavals,(boundsSing1+boundsSing2)','-o',gammavals,RLDboundSimple(gammavals)','-+')
figure
plot(gammavals,pIncCost,'-+')
