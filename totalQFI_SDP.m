function [trFIbound,hvec ] = totalQFI_SDP(KrausOps,KrausOpsDeriv,npar,nKraus,variant,dimIn,varargin)
    % solves the SDP to compute the total QFI sum_x w_X QFI_xx
    %
    % input arguments (dimensions must be self-consistent)
    % ---------------------------------------------------------------------
    % KrausOps:         all the Kraus operators stacked in columm,
    %                   expected dimensions: (dimOut*nKraus)x(dimIn)
    % KrausOpsDeriv:    the derivatives of all the Kraus operators stacked
    %                   in column as [ d1K1 ; d1K2; ..  d2K1 ; d2K2 ; .. ]
    %                   expected dimensions (dimOut*nKraus*npar)x(dimIn)
    % npar:             number of parameters
    % nKraus:           number of Kraus operators
    % version:          either 'single' (default) or 'asymptotic' 
    % dimIn:            dimension of the input Hilbert space         
    %
    % optional input arguments
    % ---------------------------------------------------------------------
    % dimOut:           dimension of the output Hilbert space, if not given
    %                   dimOut=dimIn is assumed (square Kraus operators)
    % wInv:             vector of *inverse* weights sum_x w_x QFI_xx
    %
    % output
    % ---------------------------------------------------------------------
    % trFIbound:        total QFI (scalar)
    % hvec:             optimal Hermitian matrices h_x, returned as a 
    %                   (nKraus)x(nKraus)x(npar) tensor


if nargin < 8
    wInv=ones(npar,1);
    if nargin < 7
        dimOut=dimIn;
    else
        dimOut=varargin{1};
    end
else
    dimOut=varargin{1};
    wInv=varargin{2};
end

t=sdpvar(1);
hvec=sdpvar(nKraus,nKraus,npar,'hermitian','complex');

dKtil=zeros(size(KrausOpsDeriv),'like',sdpvar);
for p=1:npar
    dKtil( 1+dimOut*nKraus*(p-1) : dimOut*nKraus*p ,:) = KrausOpsDeriv( 1+dimOut*nKraus*(p-1) : dimOut*nKraus*p ,:) - 1j*kron(hvec(:,:,p),eye(dimOut))*KrausOps;
end

Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , kron(diag(wInv),eye(dimOut*nKraus)) ] >= 0 ];
if strcmp(variant,'asymptotic')
%     eqconst = zeros( dimIn*npar,dimIn,'like',sdpvar);
    for p=1:npar
%         eqconst( 1+dimIn*(p-1) : dimIn*p , : ) = 
        Constraints = [ Constraints , KrausOps'*kron( hvec(:,:,p),speye(dimOut) )*KrausOps - 1j*KrausOpsDeriv( 1+dimOut*nKraus*(p-1) : dimOut*nKraus*p ,:)'*KrausOps == 0 ] ; 
    end
    % Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , eye(dimOut*npar*nKraus)] >= 0 , eqconst==zeros(dimIn*npar,dimIn) ];
%     Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , kron(diag(wInv),eye(dimOut*nKraus)) ] >= 0 , eqconst==zeros(dimIn*npar,dimIn) ];
elseif strcmp(variant,'single')
    % Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , eye(dimOut*npar*nKraus)] >= 0 ];
% 	Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , kron(diag(wInv),eye(dimOut*nKraus)) ] >= 0 ];
else
    error('The fifth input must be the string "single" or "asymptotic"')
end

% Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , eye(dimOut*npar*nKraus)] >= 0 ];
% Constraints = [ [ t*eye(dimIn) , dKtil' ; dKtil , kron(diag(wInv),eye(dimOut*nKraus)) ] >= 0 ];

Objective = t;
  
% options = sdpsettings('verbose',1,'solver','scs');
options = sdpsettings('verbose',1,'solver','mosek');
  
sol = optimize(Constraints,Objective,options);

trFIbound=4*value(Objective);

hvec=value(hvec);

end
