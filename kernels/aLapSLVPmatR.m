function [A,mu_sx,mu_sws] = aLapSLVPmatR(t,s)
% analytic 1-form, Laplace Single Layer Volume Potential matrix on regular box
%
% 
% Hai 06/03/21

% acquire attributes of s
p = s.p; nodeR = s.node; sc = s.sc; % order, regular nodes of a patch D, center
gridsize = abs(nodeR(1)-nodeR(2)); % regular box size...

% kernel
% logkernel = @(r0,r) -1/(4*pi)*log((r0(1)-r(1,:)).^2+(r0(2)-r(2,:)).^2); 

% some parameters 
refm = 4; reff = 2;
[x0,w0] = cheby(p); [x0m,w0m] = cheby(refm*p); [x0f,w0f] = cheby(reff*p);    
lambdac = 1.4; lambdam = 2;

% regular box [-1,1]^2
sR.Z = @(t)  sqrt(2)*(  ((t>=0)&(t<pi/2)).*(1./(sin(t)+cos(t))).*(cos(t)+1i*sin(t))*exp(1i*pi/4)...
                     + ((t>=pi/2)&(t<pi)).*(1./(sin(t-pi/2)+cos(t-pi/2))).*(cos(t-pi/2)+1i*sin(t-pi/2))*exp(1i*3*pi/4)...
                     + ((t>=pi)&(t<3*pi/2)).*(1./(sin(t-pi)+cos(t-pi))).*(cos(t-pi)+1i*sin(t-pi))*exp(1i*5*pi/4)...
                     + ((t>=3*pi/2)&(t<=2*pi)).*(1./(sin(t-3*pi/2)+cos(t-3*pi/2))).*(cos(t-3*pi/2)+1i*sin(t-3*pi/2))*exp(1i*7*pi/4));
% xToP = xToPara(sR); % maybe add this to attributes of s...

% build chebyshev approximation matrix
T = 2*cos(pi*((1:p)-1)'*(2*(p:-1:1)-1)/(2*p))/p; T(1,:) = 1/2*T(1,:); % chebyshev approximation
Tfull = zeros(p^2,p^2);
for l=1:p
    Tfull((l-1)*p+(1:p),:) = kron(eye(p),T(l,:));
end
Tfull = kron(eye(p),T)*Tfull;

% matrix from chebyshev coefficient to x^i*y*j coefficients... 
CMfull = zeros(p^2,p^2);
Ct = chebyTriCoeff(p); % relation between cheby polynomial and monomial
for k=1:p^2
    j = ceil(k/p); i = k-(j-1)*p;
    cti = Ct(i,1:i); ctj = Ct(j,1:j);
    Cm = zeros(p); Cm(1:i,1:j) = cti'*ctj;
    CMfull(:,k) = Cm(:);
end

% where mu needs to be sampled, and smooth quadr corresponding to these nodes
[muXX,muYY] = meshgrid(x0);
mu_sx = real(sc)+gridsize/2*muXX(:) + 1i*(imag(sc)+gridsize/2*muYY(:)); % where to sample density
mu_sws = w0(:)*w0(:)'; mu_sws = mu_sws(:)*(gridsize/2)^2;

% initialize matrix A
if ~isempty(t) 
    t.x = 2*(t.x-sc)/gridsize;
    A = zeros(numel(t.x),p^2); %
else % if just s, then return matrix for self and neighbor
    [tXX,tYY] = meshgrid(x0m); 
    t.x = 3*(tXX(:)+1i*tYY(:)); % 3 x 3 neighbor...
    A = zeros(numel(t.x),p^2); % 
end

% separate t.x into 3 categories
%
% close (singular quadrature)
tidxc = (abs(real(t.x))<lambdac) & (abs(imag(t.x))<lambdac);
tclo = t.x(tidxc);
t0clo = mod(mod(angle(tclo)-pi/4,2*pi)',2*pi); % 

% midway (upsampling smooth quadrture)
tidxm = (abs(real(t.x))<lambdam) & (abs(imag(t.x))<lambdam) & ~tidxc;
tmid = t.x(tidxm);

% far (smooth quadrature)
tfar = t.x(~tidxc&~tidxm);

% order of 1-form
order = 2*p-1;

%
tol = 1e-13; ratio = 1/6; refFactor = 4;
tpan0 = linspace(0,2*pi,17); [x00,w00] = gauss(ceil(3/2*p));
[YPow,XPow] = meshgrid(0:p-1); 
s0 = sR; s0.p = 2*p; s0.tpan = tpan0; s0 = quadr(s0, [], 'p', 'G'); X0 = gauss(s0.p); % a base panel for interpolation
s00 = sR; s00.p = ceil(3/2*p); s00.tpan = tpan0; s00 = quadr(s00, [], 'p', 'G'); % a reference panel for computation
r = [real(s00.x),imag(s00.x)]';
Dx0 = s00.ws.*real(s00.tang); Dy0 = s00.ws.*imag(s00.tang);
PpowerVal = (-bsxfun(@power,r(2,:),YPow(:)+1).*bsxfun(@power,r(1,:),XPow(:)));
QpowerVal = (bsxfun(@power,r(2,:),YPow(:)).*bsxfun(@power,r(1,:),XPow(:)+1));

%---------build close eval matrix-------
LapSLVPmatc = zeros(numel(tclo),p^2);
for k = 1:numel(tclo)  % loop over close targets (need to build target-dependent quadrature)
    
    r0 = [real(tclo(k));imag(tclo(k))];

    % build singular quadrature for r0
    t0 = t0clo(k);
    tlo = tpan0(sum(t0>=tpan0)); thi = tpan0(sum(t0>=tpan0)+1); 
    plen = thi - tlo;
    tmp = mod([-ratio.^(1:refFactor),ratio.^(1:refFactor)]*plen+t0,2*pi);
    tmp = tmp(tmp>0);
    if min(abs(t0-s00.tpan)) > tol % no need if tmid(j) is already in s00.tpan (base panel span)
        tpan_extra = unique([t0,tmp]);
    else
        tpan_extra = tmp;
    end
    idx_extra = sum(tpan_extra(:)'>s0.tpan(:)); uidx_extra = unique(idx_extra);
    if min(abs(uidx_extra))==0

        keyboard % why index is 0, try to get rid of this case...
    end
    idx_sx = reshape(1:numel(s00.x),s00.p,s00.np); idx_sx = idx_sx(:,setdiff(1:s00.np,uidx_extra)); idx_sx = idx_sx(:);
    s_extra.x = []; s_extra.tang = []; s_extra.ws = [];
    for j=1:numel(uidx_extra)
        tlo = s0.tpan(uidx_extra(j)); thi = s0.tpan(uidx_extra(j)+1);
        tpan_extra_k = unique([tlo,thi,tpan_extra(idx_extra==uidx_extra(j))]);
        pt = tpan_extra_k(2:end) - tpan_extra_k(1:end-1);
        t_extra0 = bsxfun(@plus,tpan_extra_k(1:end-1),bsxfun(@times,(1+x00)/2,pt)); t_extra = 2*(t_extra0(:)-(thi+tlo)/2)/(thi-tlo);
        w_extra = bsxfun(@times,w00/2,pt); 
        warning('off','MATLAB:nearlySingularMatrix');
        L = interpmat_1d(t_extra,X0);
        x_extra = L*s0.x((uidx_extra(j)-1)*s0.p+1:uidx_extra(j)*s0.p);
        xp_extra = L*s0.xp((uidx_extra(j)-1)*s0.p+1:uidx_extra(j)*s0.p); sp_extra = abs(xp_extra); 
        s_extra.tang = [s_extra.tang;xp_extra./sp_extra]; s_extra.ws = [s_extra.ws;w_extra(:).*sp_extra]; s_extra.x = [s_extra.x;x_extra]; 
    end
    sx = [s00.x(idx_sx);s_extra.x]; 
    
    % log 1-form is different for each target, therefore have to do this step repeatedly
    r = [real(sx),imag(sx)]'; 
    if order > 4, L = log1Form(r0,r,order); else, L = log1Form(r0,r,4); end
    
    r_extra = [real(s_extra.x),imag(s_extra.x)]';
    Dx = [Dx0(idx_sx);s_extra.ws.*real(s_extra.tang)]; Dy = [Dy0(idx_sx);s_extra.ws.*imag(s_extra.tang)];
    LVal = L(XPow(:)+YPow(:)+1,:);
    PVal = ([PpowerVal(:,idx_sx),-bsxfun(@power,r_extra(2,:),YPow(:)+1).*bsxfun(@power,r_extra(1,:),XPow(:))].*LVal)*Dx;
    QVal = ([QpowerVal(:,idx_sx),bsxfun(@power,r_extra(2,:),YPow(:)).*bsxfun(@power,r_extra(1,:),XPow(:)+1)].*LVal)*Dy;
    
    LapSLVPmatc(k,:) = -1/(4*pi)*(PVal' + QVal');
    
%     sR.p = ceil(3/2*p); 
%     
%     sR.tpan = [unique(mod([tpan0,t0-plen*(ratio).^(1:refFactor),t0+plen*(ratio).^(1:refFactor),t0],2*pi)),2*pi];
%     sR = quadr(sR,[],'p','G');
%     r = [real(sR.x),imag(sR.x)]';
%     Dx = sR.ws.*real(sR.tang); Dy = sR.ws.*imag(sR.tang);
%     
%     
%     % log 1-form is different for each target, therefore have to do this step repeatedly
%     if order > 4, L = log1Form(r0,r,order); else, L = log1Form(r0,r,4); end
%     
%     % assemble P*dx + Q*dy for each basis
%     LVal = L(XPow(:)+YPow(:)+1,:);
%     PpowerVal = (-bsxfun(@power,r(2,:),YPow(:)+1).*bsxfun(@power,r(1,:),XPow(:)));
%     QpowerVal = (bsxfun(@power,r(2,:),YPow(:)).*bsxfun(@power,r(1,:),XPow(:)+1));
%     PVal = (PpowerVal.*LVal)*Dx;
%     QVal = (QpowerVal.*LVal)*Dy;
%     
%     LapSLVPmatc(k,:) = -1/(4*pi)*(PVal' + QVal');
end

A(tidxc,:) = LapSLVPmatc*CMfull*Tfull;

%----------build midway eval matrix---------
u=x0; v=u;
uM = ones(p,p); for j=2:p, uM(:,j) = uM(:,j-1).*u; end
vM = ones(p,p); for j=2:p, vM(:,j) = vM(:,j-1).*v; end
Mat = kron(uM,vM);
if numel(tmid)~=0
    [sXXm,sYYm] = meshgrid(x0m); smx = sXXm(:)+1i*sYYm(:); smws = w0m(:)*w0m(:)'; smws = smws(:)';
    % build upsampling smooth quadrature for tmid
    uf = x0m; vf=uf; m=refm*p;
    ufM = ones(m,p); for j=2:p, ufM(:,j) = ufM(:,j-1).*uf; end
    vfM = ones(m,p); for j=2:p, vfM(:,j) = vfM(:,j-1).*vf; end
    Matf = kron(ufM,vfM); 
    Pm = (Mat'\Matf')';  % upsampling matrix (from mu to muf)
    LapSLVPmatm = -1/(2*pi)*bsxfun(@times,log(abs(bsxfun(@minus,tmid,smx.'))),smws); % here it is abs, therefore -1/(2*pi)
    A(tidxm,:) = LapSLVPmatm*Pm;
end

%----------build far eval matrx------------
[sXXf,sYYf] = meshgrid(x0f); sfx = sXXf(:)+1i*sYYf(:); sfws = w0f(:)*w0f(:)'; sfws = sfws(:)';
% build upsampling smooth quadrature for tfar
uf = x0f; vf=uf; m=reff*p;
ufM = ones(m,p); for j=2:p, ufM(:,j) = ufM(:,j-1).*uf; end
vfM = ones(m,p); for j=2:p, vfM(:,j) = vfM(:,j-1).*vf; end
Matf = kron(ufM,vfM); 
Pf = (Mat'\Matf')';  % upsampling matrix (from mu to muf)
if numel(tfar)~=0
    LapSLVPmatf= -1/(2*pi)*bsxfun(@times,log(abs(bsxfun(@minus,tfar,sfx.'))),sfws); % here it is abs, therefore -1/(2*pi)
    A(~tidxc&~tidxm,:) = LapSLVPmatf*Pf;
end

%----------build correction matrix for different gridsize------------
LapSLVPmatg= -1/(2*pi)*log(gridsize/2)*sfws*Pf;

%----------------finally----------------
A = (gridsize/2)^2*bsxfun(@plus,A,LapSLVPmatg);


% keyboard


end

%% 1-form
function L = log1Form(r0,r,order)
% this function returns 1-form needed up to 'order', for single target r0
% 
% it is likely for different r0, r will be different... therefore it is ok
% to call this function every time a new target comes in...
%
% however, a better way is to have an underlying fixed r for all target, do 
% adaptive quadrature, interpolation and compression based on r0...
%
%
% Hai 10/21/20

tol = 1e-13;
% recursive relation
logr = log(r(1,:).^2+r(2,:).^2);
if abs(r0(1))<tol && abs(r0(2))<tol % if the origin... much easier
    L = NaN(order,numel(r(1,:)));
    for k=1:order
        L(k,:) = 1/(k+1)*logr - 2/(k+1)^2;
    end
else    % else (there might be a complex way of setting this up, but for now, this works)
    A = abs(r0(1)*r(2,:)-r0(2)*r(1,:))./(r(1,:).^2+r(2,:).^2); B = (r0(1)*r(1,:)+r0(2)*r(2,:))./(r(1,:).^2+r(2,:).^2); % this is A = |x0y-xy0|/|r|^2; B = r0\cdot r/|r|^2
    A2 = A.^2; A2B20 = A2 + B.^2; A2B21 = A2 + (B-1).^2; logA2B21 = log(A2B21)-1; A2B21logA2B21 = A2B21.*logA2B21; % A^2; A^2+B^2; A^2+(1-B)^2; log(A^2+(1-B)^2)
    L0 = (1-B).*(logA2B21-1) + B.*(log(A2B20)-2) + 2*A.*(atan((1-B)./A)-atan(-B./A)); % 0th order
    % L02 = (1-B).*(logA2B21-1) + B.*(log(A2B20)-2) + 2*A.*atan(A./(A2B20-B)); % 0th order, there are some issue with this...
    L1 = B.*L0 + 1/2*(A2B21.*(logA2B21+1)-(1-B).^2) - 1/2*(A2B20.*log(A2B20)-B.^2); % 1st order
    L = NaN(order,numel(r(1,:))); L(1,:) = L1; L(2,:) = B.*(2*2/(2+1)*L1-2*(2-1)/(2*(2+1))) + A2B20.*(-(2-1)/(2+1)*L0+1/(2+1)) + 1/(2+1)*A2B21logA2B21 + (2-1)/(2+1)^2;
    for k=3:order
        L(k,:) = B.*(2*k/(k+1)*L(k-1,:)-2*(k-1)/(k*(k+1))) + A2B20.*(-(k-1)/(k+1)*L(k-2,:)+1/(k+1)) + 1/(k+1)*A2B21logA2B21 + (k-1)/(k+1)^2;
        L(k-2,:) = logr/(k-1)+L(k-2,:); 
    end
    L(order-1,:) = logr/order+L(order-1,:); L(order,:) = logr/(order+1)+L(order,:); 
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
