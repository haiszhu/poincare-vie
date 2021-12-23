function [A,mu_sx,mu_sws,coeffFun] = aLapSLVPmatIr(t,s)
%
%
%
% Hai 06/11/21

p = s.p;
nside = numel(s.node); nodeR = s.node0; sc = s.sc0;
tpan0 = linspace(0,2*pi,2*nside+1); s.tpan = linspace(0,2*pi,nside+1);
gridsize = abs(nodeR(1)-nodeR(2)); % regular box size...

A = zeros(numel(t.x),p^2); %
lambdac = 1.4; lambdam = 2; ratio = 1/6; refFactor = 4;
[powY,powX] = meshgrid([0:p-1]);

% find centroid...
sIr.Z = @(t) 2*(s.Z(t)-sc)/gridsize; sIr.tpan = linspace(0.0,2*pi,nside+1); sIr.p = 2*p;
[sIr, ~, ~] = quadr(sIr, [], 'p', 'G');
sIrwst = sIr.ws.*sIr.tang; sIrarea = (-1/2*(real(sIrwst))'*(imag(sIr.x)) + 1/2*(imag(sIrwst))'*(real(sIr.x)));
sIrCenx = (-1/3*(real(sIrwst))'*(real(sIr.x).*imag(sIr.x)) + 1/3*(imag(sIrwst))'*(real(sIr.x).^2))/sIrarea; sIrCeny = (-1/3*(real(sIrwst))'*(imag(sIr.x).^2) + 1/3*(imag(sIrwst))'*(imag(sIr.x).*real(sIr.x)))/sIrarea;
sIrsc = gridsize/2*(sIrCenx+1i*sIrCeny)+sc;

% rescale according to centroid
t.x = 2/gridsize*(t.x-sIrsc);
if isfield(s,'Z'),s.Z = @(t) 2/gridsize*(s.Z(t)-sIrsc); else, error('Need to provide s.Z...'); end

% ---------build approximation---------
box.Z = s.Z; box.tpan = linspace(0.0,2*pi,nside+1); box.p = p;
[box, ~, ~] = quadr(box, [], 'p', 'G'); % this is just a shift of sIr.x, one extra quadr is ok...
[x0b,w0b] = gauss(box.p); 
r = [real(box.x),imag(box.x)]';
tr_x = (x0b(:)+1)/2*(r(1,:)); tr_y = (x0b(:)+1)/2*(r(2,:));
Dx = box.ws.*real(box.tang); Dy = box.ws.*imag(box.tang);
tr_ws = 1/4*(w0b(:).*(x0b(:)+1))*(-(r(2,:)).*Dx(:)'+(r(1,:)).*Dy(:)');
x = tr_x(:)'; y = tr_y(:)'; tmpws = tr_ws(:);

% where mu needs to be sampled, and smooth quadr corresponding to these nodes
mu_sx = gridsize*x(:)/2+real(sIrsc) + 1i*(gridsize*y(:)/2+imag(sIrsc)); % where to sample density
mu_sws = tmpws*(gridsize/2)^2;

% approximation coeffcient fun
% box.p = 2*p; [box, ~, ~] = quadr(box, [], 'p', 'G'); % this is just a shift of sIr.x, one extra quadr is ok...
% [x0b,w0b] = gauss(box.p); 
% r = [real(box.x),imag(box.x)]';
% tr_x = (x0b(:)+1)/2*(r(1,:)); tr_y = (x0b(:)+1)/2*(r(2,:));
% Dx = box.ws.*real(box.tang); Dy = box.ws.*imag(box.tang);
% tr_ws = 1/4*(w0b(:).*(x0b(:)+1))*(-(r(2,:)).*Dx(:)'+(r(1,:)).*Dy(:)');
% x = tr_x(:)'; y = tr_y(:)'; tmpws = tr_ws(:);
% coeffFun =@(mu) muToCoeff(mu,x,y,tmpws,p);

N = numel(x);
Tx = NaN(p,N); Tx(1,:) = ones(1,N); Tx(2,:) = x; 
Ty = NaN(p,N); Ty(1,:) = ones(1,N); Ty(2,:) = y; 
for j=3:p
    Tx(j,:) = 2*x.*Tx(j-1,:)-Tx(j-2,:);
    Ty(j,:) = 2*y.*Ty(j-1,:)-Ty(j-2,:);
end
BLmat = NaN(p^2,p^2); %[powY,powX]=meshgrid(0:p-1); 
for i=1:p^2
    for j=i:p^2
        tmp = (Tx(powX(i)+1,:).*Ty(powY(i)+1,:).*Tx(powX(j)+1,:).*Ty(powY(j)+1,:))*tmpws;
        BLmat(i,j) = tmp; BLmat(j,i) = tmp;
    end
end

CMfull = zeros(p^2,p^2);
Ct = chebyTriCoeff(p); % relation between cheby polynomial and monomial
for l=1:p^2
    j = ceil(l/p); i = l-(j-1)*p;
    cti = Ct(i,1:i); ctj = Ct(j,1:j);
    Cm = zeros(p); Cm(1:i,1:j) = cti'*ctj;
    CMfull(:,l) = Cm(:);
end

coeffFun =@(mu) muToCoeff2(mu,tmpws,BLmat,CMfull,Tx,Ty,powX,powY,p);


tol = 1e-13; 
if isfield(s,'tcusp')
    tcusp = s.tcusp;
else
    tcusp = s.xToP(t.x); % probably want to make sure xToP is done already...
end

% distance rule (is this sufficient...)
s.node = s.Z(s.tpan(1:end-1)); 
dist = max(abs(s.node));

% separate t.x into 3 categories
% close (singular quadrature)
tidxc = abs(t.x)<lambdac*dist;
tclo = t.x(tidxc);
t0clo = tcusp(tidxc);

% midway (upsampling smooth quadrture)
tidxm = abs(t.x)<lambdam*dist;
tidxm = tidxm & ~tidxc;
tmid = t.x(tidxm);

% far (smooth quadrature)
tfar = t.x(~tidxc&~tidxm);

%---------build close eval matrix-------
if numel(tclo) > 0
    order = 2*p-1;
    % once shifted and rescaled, this irregular box is standardized to about
    % [-1,1]^2, and centered at origin
    LapSLVPmatc = zeros(numel(tclo),p^2);
    % loop over targets for singular quadr
    p0 = 3*p; [x0,w0] = gauss(p0);
    s.tpan = tpan0; 
    s0 = s; s0.p = ceil(3/2*p0); s0 = quadr(s0, [], 'p', 'G'); X0 = gauss(s0.p); % a base panel for interpolation
    s00 = s; s00.p = p0; s00 = quadr(s00, [], 'p', 'G'); % a reference panel for computation
    % figure out all extra nodes and their corresponding index to target tclo(j)
    Xt_extra = cell(1,numel(tclo)); Xw_extra = cell(1,numel(tclo)); 
    Idx_sx = cell(1,numel(tclo)); Idx_sx_extra = cell(1,numel(tclo)); 
    idx_start = 0;
    for j = 1:numel(tclo)
        % refFactor = 4; ratio = 1/6;
        tmp = mod([-ratio.^(1:refFactor),ratio.^(1:refFactor)]*2*pi/numel(tpan0)+t0clo(j),2*pi);
        tmp = tmp(tmp>0); % exclude 0...
        if min(abs(t0clo(j)-s00.tpan)) > tol % no need if tmid(j) is already in s00.tpan (base panel span)
            tpan_extra = unique([t0clo(j),tmp]);
        else
            tpan_extra = tmp;
        end
        idx_extra = sum(tpan_extra(:)'>s0.tpan(:)); uidx_extra = unique(idx_extra);
        idx_sx = reshape(1:numel(s00.x),s00.p,s00.np); idx_sx = idx_sx(:,setdiff(1:s00.np,uidx_extra)); idx_sx = idx_sx(:);
        xt_extra = []; xw_extra = [];
        for k=1:numel(uidx_extra)
            tlo = s0.tpan(uidx_extra(k)); thi = s0.tpan(uidx_extra(k)+1);
            tpan_extra_k = unique([tlo,thi,tpan_extra(idx_extra==uidx_extra(k))]);
            pt = tpan_extra_k(2:end) - tpan_extra_k(1:end-1);
            t_extra0 = bsxfun(@plus,tpan_extra_k(1:end-1),bsxfun(@times,(1+x0)/2,pt)); 
            w_extra = bsxfun(@times,w0/2,pt); 
            xt_extra = [xt_extra;t_extra0(:)]; xw_extra = [xw_extra;w_extra(:)];
        end
        Xt_extra{j} = xt_extra; Xw_extra{j} = xw_extra; 
        tmp = numel(xt_extra);
        Idx_sx_extra{j} =  idx_start+(1:tmp);
        Idx_sx{j} = idx_sx;
        idx_start = idx_start + tmp;
    end
    % tangent and quadrature weight at these extra nodes
    Xt_extra = vertcat(Xt_extra{:}); Xw_extra = vertcat(Xw_extra{:}); 
    X_extra = s.Z(Xt_extra); Xp_extra = zeros(size(X_extra));
    for k=1:numel(s0.tpan)-1
        tlo = s0.tpan(k); thi = s0.tpan(k+1);
        t0_idx = (Xt_extra>=tlo)&(Xt_extra<=thi);
        t0_extra = Xt_extra(t0_idx);
        warning('off','MATLAB:nearlySingularMatrix');
        Lmat = interpmat_1d(2*(t0_extra-(tlo+thi)/2)/(thi-tlo),X0);
        warning('on','MATLAB:nearlySingularMatrix')
        Xp_extra(t0_idx) = Lmat*s0.xp((k-1)*s0.p+(1:s0.p));
    end
    Xtang_extra = Xp_extra./abs(Xp_extra); Xws_extra = Xw_extra.*abs(Xp_extra);
    % 1-form, each target has one matrix of dimension (order=2*p-1) X (number of corresponding source points along s.Z)
    Ltclo = cell(1,numel(tclo));
    if order <= 4, order = 4; end
    for j = 1:numel(tclo)
        Ltclo{j} = log1Form([real(tclo(j));imag(tclo(j))],[[real(s00.x(Idx_sx{j}));real(X_extra(Idx_sx_extra{j}))],[imag(s00.x(Idx_sx{j}));imag(X_extra(Idx_sx_extra{j}))]]',order);
    end
    % actually singular quadrature computation (after all previous preparation...) 
    r = [real(s00.x),imag(s00.x)]';
    Dx0 = s00.ws.*real(s00.tang); Dy0 = s00.ws.*imag(s00.tang);
    PpowerVal0 = (-bsxfun(@power,r(2,:),powY(:)+1).*bsxfun(@power,r(1,:),powX(:)));
    QpowerVal0 = (bsxfun(@power,r(2,:),powY(:)).*bsxfun(@power,r(1,:),powX(:)+1));
    r_extra = [real(X_extra),imag(X_extra)]';
    Dx_extra = Xws_extra.*real(Xtang_extra); Dy_extra = Xws_extra.*imag(Xtang_extra);
    XpowerVal_extra = bsxfun(@power,r_extra(1,:),(0:p)'); YpowerVal_extra = bsxfun(@power,r_extra(2,:),(0:p)'); 
    PpowerVal_extra = -YpowerVal_extra(powY(:)+2,:).*XpowerVal_extra(powX(:)+1,:);
    QpowerVal_extra =  YpowerVal_extra(powY(:)+1,:).*XpowerVal_extra(powX(:)+2,:);
    for j = 1:numel(tclo)
        % log 1-form is different for each target, therefore have to do this step repeatedly
        LVal = Ltclo{j}(powX(:)+powY(:)+1,:);
        Dx = [Dx0(Idx_sx{j});Dx_extra(Idx_sx_extra{j})]; Dy = [Dy0(Idx_sx{j});Dy_extra(Idx_sx_extra{j})];
        PVal = ([PpowerVal0(:,Idx_sx{j}),PpowerVal_extra(:,Idx_sx_extra{j})].*LVal)*Dx;
        QVal = ([QpowerVal0(:,Idx_sx{j}),QpowerVal_extra(:,Idx_sx_extra{j})].*LVal)*Dy;
        LapSLVPmatc(j,:) = -1/(4*pi)*(PVal(:) + QVal(:))';
    end
    
    % smooth quadr for scaling
    s.tpan = linspace(0,2*pi,nside+1); s.p = p; 
    [s, ~, ~] = quadr(s, [], 'p', 'G'); swst = s.ws.*s.tang; 
    r = [real(s.x),imag(s.x)]'; 
    [t0c,w0c] = gauss(p); t0c = (t0c+1)/2; w0c = w0c/2;
    tr_x = t0c*r(1,:); tr_y = t0c*r(2,:);  % t\in [0,1] for 1-form construction
    Dx = real(swst); Dy = imag(swst);
    tr_ws = (w0c(:).*t0c(:))*(-(r(2,:)).*Dx(:)'+(r(1,:)).*Dy(:)');
    XpowerVal_smooth = bsxfun(@power,tr_x(:)',(0:p-1)'); YpowerVal_smooth = bsxfun(@power,tr_y(:)',(0:p-1)'); 
    muTxy = XpowerVal_smooth(powX(:)'+1,:).*YpowerVal_smooth(powY(:)'+1,:);
    LapSLVPmatc = (gridsize/2)^2*bsxfun(@plus,LapSLVPmatc,- 1/(2*pi)*log(gridsize/2)*(muTxy*tr_ws(:))');

    A(tidxc,:) = LapSLVPmatc;
end
    
%----------build midway eval matrix--------
if numel(tmid)>0
    s.tpan = linspace(0,2*pi,nside+1); s.p = 3*p; 
    [sm, ~, ~] = quadr(s, [], 'p', 'G'); swst = sm.ws.*sm.tang; 
    rm = [real(sm.x),imag(sm.x)]'; 
    [t0m,w0m] = gauss(3*p); t0m = (t0m+1)/2; w0m = w0m/2;
    tr_x = t0m*rm(1,:); tr_y = t0m*rm(2,:);  % t\in [0,1] for 1-form construction
    Dx = real(swst); Dy = imag(swst);
    tr_ws = (w0m(:).*t0m(:))*(-(rm(2,:)).*Dx(:)'+(rm(1,:)).*Dy(:)');
    smx = tr_x(:)'+1i*tr_y(:)'; smws = tr_ws(:)'; 
    
    XpowerValm = bsxfun(@power,real(smx),(0:p-1)'); YpowerValm = bsxfun(@power,imag(smx),(0:p-1)'); 
    muTxy = XpowerValm(powX(:)'+1,:).*YpowerValm(powY(:)'+1,:);
    LapSLVPmatm = - 1/(2*pi)*log(abs(bsxfun(@minus,tmid(:),smx)))*bsxfun(@times,muTxy',smws(:));
    LapSLVPmatm = (gridsize/2)^2*bsxfun(@plus,LapSLVPmatm,- 1/(2*pi)*log(gridsize/2)*(muTxy*smws(:))');
    A(tidxm,:) = LapSLVPmatm;
end

%----------build far eval matrx------------
if numel(tfar)>0
    s.tpan = linspace(0,2*pi,nside+1); s.p = 2*p; 
    [sf, ~, ~] = quadr(s, [], 'p', 'G'); swst = sf.ws.*sf.tang; 
    rf = [real(sf.x),imag(sf.x)]'; 
    [t0f,w0f] = gauss(2*p); t0f = (t0f+1)/2; w0f = w0f/2;
    tr_x = t0f*rf(1,:); tr_y = t0f*rf(2,:);  % t\in [0,1] for 1-form construction
    Dx = real(swst); Dy = imag(swst);
    tr_ws = (w0f(:).*t0f(:))*(-(rf(2,:)).*Dx(:)'+(rf(1,:)).*Dy(:)');
    sfx = tr_x(:)'+1i*tr_y(:)'; sfws = tr_ws(:)'; 
    
    XpowerValf = bsxfun(@power,real(sfx),(0:p-1)'); YpowerValf = bsxfun(@power,imag(sfx),(0:p-1)'); 
    muTxy = XpowerValf(powX(:)'+1,:).*YpowerValf(powY(:)'+1,:);
    LapSLVPmatf = - 1/(2*pi)*log(abs(bsxfun(@minus,tfar(:),sfx)))*bsxfun(@times,muTxy',sfws(:));
    LapSLVPmatf = (gridsize/2)^2*bsxfun(@plus,LapSLVPmatf,- 1/(2*pi)*log(gridsize/2)*(muTxy*sfws(:))');
    A(~tidxc&~tidxm,:) = LapSLVPmatf;
end


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

function coeff = muToCoeff2(mutxy,tmpws,BLmat,CMfull,Tx,Ty,powX,powY,p)
% from mu to polynomial approximation coefficients...
%
%
% Hai 06/12/21

% Lmat = interpmat_1d(gauss(2*p),gauss(p));
% mutmp = reshape(mutxy,p,numel(mutxy)/p);
% mutxy1 = Lmat*mutmp; mutxyf = zeros(2*p,2*numel(mutxy)/p);
% for k=1:(numel(mutxy)/p^2)
%     mutxyf(:,(k-1)*2*p+(1:2*p)) = mutxy1(:,(k-1)*p+(1:p))*Lmat';
% end
% xtmp = reshape(x,2*p,numel(x)/(2*p));
% ytmp = reshape(y,2*p,numel(y)/(2*p));

BLrhs = NaN(p^2,1); mutmp = mutxy(:);
for i=1:p^2 
    BLrhs(i) = (Tx(powX(i)+1,:).*Ty(powY(i)+1,:))*(mutmp.*tmpws); 
end
warning('off','MATLAB:nearlySingularMatrix');
C_cheb = BLmat\BLrhs; C_cheb(abs(C_cheb)<=1e-09) = 0;
warning('on','MATLAB:nearlySingularMatrix')
coeff = CMfull*C_cheb;

end

%% approximation coefficient
function coeff = muToCoeff(mutxy,x,y,tmpws,p)
% from mu to polynomial approximation coefficients...
%
%
% Hai 06/12/21

% Lmat = interpmat_1d(gauss(2*p),gauss(p));
% mutmp = reshape(mutxy,p,numel(mutxy)/p);
% mutxy1 = Lmat*mutmp; mutxyf = zeros(2*p,2*numel(mutxy)/p);
% for k=1:(numel(mutxy)/p^2)
%     mutxyf(:,(k-1)*2*p+(1:2*p)) = mutxy1(:,(k-1)*p+(1:p))*Lmat';
% end
% xtmp = reshape(x,2*p,numel(x)/(2*p));
% ytmp = reshape(y,2*p,numel(y)/(2*p));

N = numel(x);
Tx = NaN(p,N); Tx(1,:) = ones(1,N); Tx(2,:) = x; 
Ty = NaN(p,N); Ty(1,:) = ones(1,N); Ty(2,:) = y; 
for j=3:p
    Tx(j,:) = 2*x.*Tx(j-1,:)-Tx(j-2,:);
    Ty(j,:) = 2*y.*Ty(j-1,:)-Ty(j-2,:);
end
BLmat = NaN(p^2,p^2); [powY,powX]=meshgrid(0:p-1); 
for i=1:p^2
    for j=i:p^2
        tmp = (Tx(powX(i)+1,:).*Ty(powY(i)+1,:).*Tx(powX(j)+1,:).*Ty(powY(j)+1,:))*tmpws;
        BLmat(i,j) = tmp; BLmat(j,i) = tmp;
    end
end
BLrhs = NaN(p^2,1); mutmp = mutxy(:);
for i=1:p^2 
    BLrhs(i) = (Tx(powX(i)+1,:).*Ty(powY(i)+1,:))*(mutmp.*tmpws); 
end
warning('off','MATLAB:nearlySingularMatrix');
C_cheb = BLmat\BLrhs; C_cheb(abs(C_cheb)<=1e-09) = 0;
warning('on','MATLAB:nearlySingularMatrix')
CMfull = zeros(p^2,p^2);
Ct = chebyTriCoeff(p); % relation between cheby polynomial and monomial
for l=1:p^2
    j = ceil(l/p); i = l-(j-1)*p;
    cti = Ct(i,1:i); ctj = Ct(j,1:j);
    Cm = zeros(p); Cm(1:i,1:j) = cti'*ctj;
    CMfull(:,l) = Cm(:);
end
coeff = CMfull*C_cheb;

end


%% interpolation matrix
function L = interpmat_1d(t,s)
% INTERPMAT_1D   interpolation matrix from nodes in 1D to target nodes
%
% L = interpmat_1d(t,s) returns interpolation matrix taking values on nodes s
%  to target nodes t. Computed in Helsing style.

% bits taken from qplp/interpmatrix.m from 2013
% Barnett 7/17/16

if nargin==0, test_interpmat_1d; return; end

p = numel(s); q = numel(t); s = s(:); t = t(:);       % all col vecs
n = p; % set the polynomial order we go up to (todo: check why bad if not p)
V = ones(p,n); for j=2:n, V(:,j) = V(:,j-1).*s; end   % polyval matrix on nodes
R = ones(q,n); for j=2:n, R(:,j) = R(:,j-1).*t; end   % polyval matrix on targs
L = (V'\R')'; % backwards-stable way to do it (Helsing) See corners/interpdemo.m
end 

%% chebyshev polynomial coefficient in terms of monomial
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



