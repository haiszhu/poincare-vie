function [Cm,xyPow] = monoApprox(Fxy)
% this function is built on Chebyshev approximation..., and returns a
% monomial approximation coefficients and their corresponding power... 
% 
% it requires: cheby.m, chebyCoeff2d.m, chebyTriCoeff.m, monoCoeff2d.m
%
% it also requires Fxy to be order this way: 
%           p = 16; x = cheby(p); [XX,YY] = meshgrid(x); 
%           Fxy = reshape(fxy(XX(:)',YY(:)'),p,p);
%
% it is an approximation on regular box, therefore order is implied...
%
% it cuts insignificant coefficients smaller than 1e-14...
%
% Hai 10/21/20

N = sqrt(numel(Fxy)); % implied order
[powY,powX]=meshgrid(0:N-1);
Cmfull = monoCoeff2d(Fxy); idx = find(abs(Cmfull)>1e-14); 
Cm = Cmfull(idx)'; xPow = powX(idx); yPow = powY(idx); xyPow = [xPow';yPow'];
% fxy \approx @(x,y) Cm*(bsxfun(@power,x,xyPow(1,:)').*bsxfun(@power,y,xyPow(2,:)'));


end

function Cm = monoCoeff2d(Fxy)
% for the purpose of computing monomial interpolation coefficients
%
%   [XX,YY] = meshgrid(x);  
%   Fxy = fxy(XX,YY);  
%   fxy \approx \sum_{i,j} C(i,j)*t^(i-1)*s^(j-1);
%     
%  ** need this for potential 1-form construction
%
% Hai 09/21/20

N = sqrt(numel(Fxy)); % implied order
Ctmp = chebyCoeff2d(Fxy); % chebyshev approximation
C = zeros(size(Ctmp)); idx = find(abs(Ctmp)>1e-14); C(idx) = Ctmp(idx);
Ct = chebyTriCoeff(N); % relation between cheby polynomial and monomial
Cm = zeros(N);
for i=1:N
    for j=1:N
        cti = Ct(i,1:i); ctj = Ct(j,1:j);
        Cij = C(i,j)*cti'*ctj;
        Cm(1:i,1:j) = Cm(1:i,1:j) + Cij;
    end
end

end

function [C,T] = chebyCoeff2d(Fxy)
% for the purpose of computing chebyshev interpolation coefficients
% get chebyshev polynomial coefficients
%
% ** constant term also rescaled...
%
%   [XX,YY] = meshgrid(x);  
%   Fxy = fxy(XX,YY);  
%   fxy \approx \sum_{i,j} C(i,j)*cheb_i(t)*cheb_j(s);
%     
%
% Hai 09/21/20, modifided 05/29/21 to include the full matrix to compute
% approximation coefficient "kron(eye(N),T)*Tfull"

N = sqrt(numel(Fxy)); % asssuming same order along both direction

% form coefficient computation matrix by Clenshaw (Numerical Methods for Special Functions)
T = 2*cos(pi*((1:N)-1)'*(2*(N:-1:1)-1)/(2*N))/N; T(1,:) = 1/2*T(1,:);
C = T*Fxy; % this is equivalent to reshape(kron(eye(N),T)*Fxy(:),N,N);
C = T*C'; % C = C';
idx = abs(C)<=1e-14; C(idx) = 0;
if 0
    Tfull = zeros(N^2,N^2);
    for k=1:N
        Tfull((k-1)*N+(1:N),:) = kron(eye(N),T(k,:));
    end
    C = (kron(eye(N),T)*Tfull)*Fxy(:);
end

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
