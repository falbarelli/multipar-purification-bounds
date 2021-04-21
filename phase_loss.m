eta=.329;
%%% the analytical results, for comparison
lossBound=1/(eta*(1-eta));
phaseBoundSing=4.*eta/((1+sqrt(eta)).^2);
phaseBoundAsymp=(4*eta/(1-eta));
incompBound=2*((1-eta)./( eta + sqrt(eta) - sqrt(2*(1+sqrt(eta))))).^2;

nKraus=2;
dimIn=2;
dimOut=3;
K0= [ sqrt(eta) , 0 ; 0 , 1 ; 0 , 0];
K1= [ 0 , 0 ; 0, 0 ; sqrt(1-eta), 0 ];
KrausOps = [ K0 ; K1];

detaK0=0.5*[ 1/sqrt(eta) , 0 ; 0 , 0  ; 0 , 0];
detaK1=-0.5*[ 0 , 0 ; 0, 0 ; 1/sqrt(1-eta) , 0 ];
detaKrausOps=[ detaK0 ; detaK1 ];

genPhi=[ 1 , 0 ; 0 , 0]; 
dPhiKrausOps= -1j*[ K0*genPhi ; K1*genPhi ];

%%% single-parameter results

[ lossBoundAsymNum , hlossAs ]    = totalQFI_SDP(KrausOps,detaKrausOps,1,nKraus,'asymptotic',dimIn,dimOut);
[ lossBoundSingNum , hlossSing ]  = totalQFI_SDP(KrausOps,detaKrausOps,1,nKraus,'single',dimIn,dimOut);

[ phaseBoundAsymNum , hphaseAs ]   =totalQFI_SDP(KrausOps,dPhiKrausOps,1,nKraus,'asymptotic',dimIn,dimOut);
[ phaseBoundSingNum , hphaseSing ] =totalQFI_SDP(KrausOps,dPhiKrausOps,1,nKraus,'single',dimIn,dimOut);

%%% check that they correspond to the analytical results, this should all
%%% be zero
lossBoundAsymNum-lossBound
lossBoundSingNum-lossBound
phaseBoundAsymNum-phaseBoundAsymp
phaseBoundSingNum-phaseBoundSing

%%% 2-parameter results

KrausOpsDeriv2par=[ dPhiKrausOps ; detaKrausOps ] ;

incompSing = 2/totalQFI_SDP(KrausOps,KrausOpsDeriv2par,2,nKraus,'single',dimIn,dimOut,[phaseBoundSingNum,lossBoundSingNum]);

%%% check the analytical results, this is zero
incompBound - incompSing
%%% and this is always 1
[incompAs , hAsymp] = totalQFI_SDP(KrausOps,KrausOpsDeriv2par,2,nKraus,'asymptotic',dimIn,dimOut,[phaseBoundAsymNum,lossBoundAsymNum]);
2/incompAs

[boundAs , hAsymp2] = totalQFI_SDP(KrausOps,KrausOpsDeriv2par,2,nKraus,'asymptotic',dimIn,dimOut);
