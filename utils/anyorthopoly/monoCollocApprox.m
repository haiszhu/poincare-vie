function [Cm,xyPow] = monoCollocApprox(tr,tws,mutr,p,lambda,sCen)
% collocation approximation...
% 
% tr:   t*r, t \in [0,1]
% tws:  quadrature weight for smooth volume integral
% mutr: body force value at t*r...
% p:    collocation approximation order...
%
%
% Hai 01/11/21

if nargin < 5, lambda = 1; end
if nargin < 6, sCen = [0;0]; end

tr_x = 1/lambda*(tr(1,:)-sCen(1)); tr_y = 1/lambda*(tr(2,:)-sCen(2));
ChebTr_x = ones(p,numel(tr_x)); ChebTr_y = ChebTr_x; ChebTr_x(2,:) = tr_x; ChebTr_y(2,:) = tr_y;
for k=3:p
    ChebTr_x(k,:) = 2*tr_x.*ChebTr_x(k-1,:)-ChebTr_x(k-2,:);
    ChebTr_y(k,:) = 2*tr_y.*ChebTr_y(k-1,:)-ChebTr_y(k-2,:);
end

% collocation matrix A
A = zeros(p^2,p^2); [powY,powX]=meshgrid(0:p-1);
for i=1:p^2
    A(i,:) = 1/lambda^2*((bsxfun(@times,ChebTr_x(powX(i)+1,:).*ChebTr_y(powY(i)+1,:),ChebTr_x(powX(:)+1,:).*ChebTr_y(powY(:)+1,:)))*tws(:))';
end

% rhs
muws = 1/lambda^2*mutr(:).*tws(:);
rhs = (ChebTr_x(powX(:)+1,:).*ChebTr_y(powY(:)+1,:)*muws(:));

%
coeff = A\rhs; tmp = reshape(coeff,p,p);
C = zeros(size(tmp)); idx = find(abs(tmp)>1e-09); C(idx) = tmp(idx);
Cmfull = zeros(p); Ct = chebyTriCoeff(p); 
for i=1:p
    for j=1:p
        cti = Ct(i,1:i); ctj = Ct(j,1:j);
        Cij = C(i,j)*cti'*ctj;
        Cmfull(1:i,1:j) = Cmfull(1:i,1:j) + Cij;
    end
end

idx = find(abs(Cmfull)>1e-09); 
Cm = Cmfull(idx)'; xPow = powX(idx); yPow = powY(idx); xyPow = [xPow';yPow'];

end

function Ct = chebyTriCoeff(n)
%
%
% Hai 09/21/20 modified on 01/13/21 for order n higher than 16...

if nargin==0, test_chebyTriCoeff; return; end

if n<= 16
    Ct = zeros(16);
    Ct(1,1)       =   1;
    Ct(2,2)       =   1;
    Ct(3,1:2:3)   = [-1,2];
    Ct(4,2:2:4)   = [-3,4];
    Ct(5,1:2:5)   = [1,-8,8];
    Ct(6,2:2:6)   = [5,-20,16];
    Ct(7,1:2:7)   = [-1,18,-48,32];
    Ct(8,2:2:8)   = [-7,56,-112,64];
    Ct(9,1:2:9)   = [1,-32,160,-256,128];
    Ct(10,2:2:10) = [9,-120,432,-576,256];
    Ct(11,1:2:11) = [-1,50,-400,1120,-1280,512];
    Ct(12,2:2:12) = [-11,220,-1232,2816,-2816,1024];
    Ct(13,1:2:13) = [1,-72,840,-3584,6912,-6144,2048];
    Ct(14,2:2:14) = [13,-364,2912,-9984,16640,-13312,4096];
    Ct(15,1:2:15) = [-1,98,-1568,9408,-26880,39424,-28672,8192];
    Ct(16,2:2:16) = [-15,560,-6048,28800,-70400,92160,-61440,16384];
    Ct = Ct(1:n,1:n);
else
    Ct = zeros(n);
    Ct(1,1)       =   1;
    Ct(2,2)       =   1;
    for k=3:n
        Ct(k,1:k) = [0,2*Ct(k-1,1:k-1)]-[Ct(k-2,1:k-2),0,0];
    end
%     keyboard
    % define for higher order... check if this is correct with Matlab built
    % in cheby fun
end

end
