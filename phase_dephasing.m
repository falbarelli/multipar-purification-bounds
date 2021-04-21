eta=.438;
sx=full(Pauli(1));
sy=full(Pauli(2));
sz=full(Pauli(3));

etaBound=1/(1-eta^2); % noise estimation, single and asymptotic result
phiBoundAsymp=eta^2/(1-eta^2); % phase estimation, asymptotic result
phiBoundSing=eta^2; % phase estimation, single-shot result

hphiSingAnal=sx*sqrt(1-eta^2)/2+eye(2)/2
hphiAsympAnal=sx*1/(2*sqrt(1-eta^2))+eye(2)/2

nKraus=2;
dimIn=2;
dimOut=dimIn;
npar=2;

%%% single-parameter results 

K0=sqrt((1+eta)/2)*eye(2);
K1=sqrt((1-eta)/2)*sz;
KrausOps=[ K0 ; K1 ];

dphiK0=K0*diag([1j,0]);
dphiK1=K1*diag([1j,0]);
% dphiK0=K0*1j*sz/2;
% dphiK1=K1*1j*sz/2;
KrausOpsDerivPhi=[ dphiK0 ; dphiK1 ];

[ phiBoundSingNum , hphiSing ] = totalQFI_SDP(KrausOps,KrausOpsDerivPhi,1,nKraus,'single',dimIn,dimOut)
[ phiBoundAsympNum , hphiAsymp ] = totalQFI_SDP(KrausOps,KrausOpsDerivPhi,1,nKraus,'asymptotic',dimIn,dimOut)

% check of the analytical results
phiBoundSingNum-phiBoundSing
hphiSing-hphiSingAnal
phiBoundAsympNum-phiBoundAsymp
hphiAsymp-hphiAsympAnal

detaK0=eye(2)/(2*sqrt(2)*sqrt(1+eta));
detaK1=-sz/(2*sqrt(2)*sqrt(1-eta));
KrausOpsDerivEta=[ detaK0 ; detaK1 ];

[ etaBoundSingNum , hetaSing ] = totalQFI_SDP(KrausOps,KrausOpsDerivEta,1,nKraus,'single',dimIn,dimOut)
[ etaBoundAsympNum , hetaAsymp ] = totalQFI_SDP(KrausOps,KrausOpsDerivEta,1,nKraus,'asymptotic',dimIn,dimOut)

etaBoundSingNum-etaBound
etaBoundAsympNum-etaBound

%%% 2-parameter results 
KrausOpsDeriv2par=[KrausOpsDerivPhi ; KrausOpsDerivEta] ;

[ trfi2parSing , hmatsSing ] = totalQFI_SDP(KrausOps,KrausOpsDeriv2par,npar,nKraus,'single',dimIn,dimOut)
[ trfi2parAsymp , hmatsAsymp ] = totalQFI_SDP(KrausOps,KrausOpsDeriv2par,npar,nKraus,'asymptotic',dimIn,dimOut)

%%% check, these should be zeros
hmatsSing(:,:,1)-hphiSingAnal
hmatsAsymp(:,:,1)-hphiAsympAnal

trfi2parSing-(etaBound+phiBoundSing)
trfi2parAsymp-(etaBound+phiBoundAsymp)

[ incompatSing , hmatsSing ] = totalQFI_SDP(KrausOps,KrausOpsDeriv2par,npar,nKraus,'single',dimIn,dimOut,[phiBoundSingNum,etaBoundSingNum]);
[ incompatAsymp , hmatsAsymp ] = totalQFI_SDP(KrausOps,KrausOpsDeriv2par,npar,nKraus,'asymptotic',dimIn,dimOut,[phiBoundAsympNum,etaBoundAsympNum]);

%%% these should always be 1
2/incompatSing
2/incompatAsymp