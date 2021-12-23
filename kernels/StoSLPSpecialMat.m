function A = StoSLPSpecialMat(t,s,side)
% special quadrature for Stokes BVP...
% 
%
% Hai 02/28/21, modified from the implementation in "Solution of Stokes flow in complex nonsmooth 2D geometries via a linear-scaling high-order adaptive integral equation scheme "

if nargin<3, side = 'i'; end

be = 2; qntype='G'; 
Ag = StoSLPmat(t,s); A = Ag(1:end/2,:) + 1i*Ag(end/2+1:end,:);
ls=0; len=s.p; np = s.np;
panlen = zeros(1,np);
Imn = interpmat(s.p,be*s.p,qntype);
for k=1:np
    ss.p=s.p; ss.x=s.x(ls+1:ls+len); ss.t=s.t(ls+1:ls+len);
    ss.xp=s.xp(ls+1:ls+len); ss.xpp=s.xpp(ls+1:ls+len); ss.sp=s.sp(ls+1:ls+len);
    ss.w=s.w(ls+1:ls+len); ss.nx=s.nx(ls+1:ls+len);
    ss.ws=s.ws(ls+1:ls+len); ss.wxp=s.wxp(ls+1:ls+len);
    ss.tlo=s.tlo(k); ss.thi=s.thi(k);
    
    j = (k-1)*s.p+(1:s.p); 
    panlen(k) = sum(s.ws(j)); 
    ik = (abs(t.x - s.xlo(k)) + abs(t.x - s.xhi(k))) < 1.7*panlen(k);
    if sum(ik)>0
        tt.x=t.x(ik(:)); 
        ssf = quadr_panf(ss,be,qntype); ssf.t = gauss(be*ss.p);  % upsampled panel
        [As, A1, A2] = Sspecialquad(tt,ssf,s.xlo(k),s.xhi(k),side);
        S11p1iS21 = (As/2+A1.*bsxfun(@plus,-tt.x,ssf.x.')/2)*Imn;
        S12p1iS22 = (1i*As/2+A2.*bsxfun(@plus,-tt.x,ssf.x.')/2)*Imn;
        Acorr = [S11p1iS21,S12p1iS22];
        A(ik,[j,end/2+j]) = Acorr;
    end
    ls=ls+len;
end
A = [real(A);imag(A)];

end

function sf = quadr_panf(s, be, qntype)  
% set up quadrature on a closed segment
% QUADR_panf - set up quadrature (either coarse or fine nodes) on a segment struct
%
% sf = quadr_panf(s, be) gives quadrature on coarse or fine nodes.
% Inputs: s  = segment struct containing parametrization
%         be = factor by which to increase panel nodes
%  
% Outputs: sf - segment struct on fine node
if be == 1
    sf = s;
    if ~isfield(s,'p')
    s.p=16; 
    end
    p = s.p; % default panel order
    if qntype=='G', [~, w, D] = gauss(p); else, [~, w, D] = cheby(p); end 
    sf.xp = D*sf.x;
    sf.xpp = D*sf.xp;   % acceleration Z''(sf.x)
    sf.w = w;
    sf.sp = abs(sf.xp); sf.tang = sf.xp./sf.sp; sf.nx = -1i*sf.tang;    % outward unit norma
else
if ~isfield(s,'p')
    s.p=16; 
end
p = s.p; % default panel order
sf=[]; sf.p=ceil(be*s.p); pf=sf.p;
Imn = interpmat(p, pf, qntype);
sf.x = Imn*s.x;
if qntype=='G', [xx, w, D] = gauss(pf); else, [xx, w, D] = cheby(pf); end 
if ~isfield(s,'Zp') 
    if ~isfield(s,'xp'), sf.xp = D*sf.x;  else, sf.xp = Imn*s.xp*(s.thi-s.tlo)/2; end  % velocities Z'(sf.x)
else
    sf.xp = 1/2*(s.thi-s.tlo)*s.Zp(s.tlo + (1+xx)/2*(s.thi-s.tlo));
end
if ~isfield(s,'Zpp')
    if ~isfield(s,'xpp'), sf.xpp = D*sf.xp;  else, sf.xpp = Imn*s.xpp*(s.thi-s.tlo)/2; end  % acceleration Z''(sf.x)
else
    sf.xpp = 1/2*(s.thi-s.tlo)*s.Zpp(s.tlo + (1+xx)/2*(s.thi-s.tlo));
end
sf.w = w;
sf.sp = abs(sf.xp); sf.tang = sf.xp./sf.sp; sf.nx = -1i*sf.tang;    % outward unit normals
sf.cur = -real(conj(sf.xpp).*sf.nx)./sf.sp.^2;
sf.ws = sf.w.*sf.sp; % speed weights
sf.wxp = sf.w.*sf.xp; % complex speed weights (Helsing's wzp)
end
end



function [A, A1, A2] = Sspecialquad(t,s,a,b,side)
% SSPECIALQUAD - SLP val+grad close-eval Helsing "special quadrature" matrix
%
% [A] = Sspecialquad(t,s,a,b) returns
% [A An] = Sspecialquad(t,s,a,b) also gives target normal-derivs (needs t.nx)
% [A A1 A2] = Sspecialquad(t,s,a,b) also gives target x,y-derivs
% Inputs: t = target seg struct (with column-vec t.x targets in complex plane)
%         s = src node seg struct (with s.x, s.wxp, s.nx)
%         a = panel start, b = panel end, in complex plane.
% Output: A (n_targ * n_src) is source-to-target value matrix
%         An or A1, A2 = source to target normal-deriv (or x,y-deriv) matrices
% Efficient only if multiple targs, since O(p^3).
% See Helsing-Ojala 2008 (special quadr Sec 5.1-2), Helsing 2009 mixed (p=16),
% and Helsing's tutorial demo11b.m LGIcompRecFAS()
if nargin<5, side = 'i'; end     % interior or exterior
zsc = (b-a)/2; zmid = (b+a)/2; % rescaling factor and midpoint of src segment
y = (s.x-zmid)/zsc; x = (t.x-zmid)/zsc;  % transformed src nodes, targ pts
N = numel(x);                            % # of targets
p = numel(s.x);                          % assume panel order is # nodes
c = (1-(-1).^(1:p))./(1:p);              % Helsing c_k, k = 1..p.
V = ones(p,p); for k=2:p, V(:,k) = V(:,k-1).*y; end  % Vandermonde mat @ nodes
P = zeros(p+1,N);      % Build P, Helsing's p_k vectorized on all targs...
d = 1.1; inr = abs(x)<=d; ifr = abs(x)>d;      % near & far treat separately
% compute P up to p+1 instead of p as in DLP, since q_k needs them:
%gam = 1i;  % original choice: branch cut is semicircle behind panel
gam = exp(1i*pi/4);  % smaller makes cut closer to panel. barnett 4/17/18
% note gam = 1 fails, and gam = -1 put cuts in the domain.
if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
P(1,:) = log(gam) + log((1-x)./(gam*(-1-x)));  % init p_1 for all targs int

% upwards recurrence for near targets, faster + more acc than quadr...
% (note rotation of cut in log to -Im; so cut in x space is lower unit circle)
Nn =  numel(find(inr));
if Nn ~= 0 
    for k=1:p
        P(k+1,inr) = x(inr).'.*P(k,inr) + c(k); 
    end  % recursion for p_k
end
% fine quadr (no recurrence) for far targets (too inaccurate cf downwards)...
Nf =  numel(find(ifr)); wxp = s.wxp/zsc; % rescaled complex speed weights
if Nf>0 
    P(end,ifr) = sum(((wxp.*(V(:,end).*y(:)))*ones(1,Nf))./bsxfun(@minus,y,x(ifr).'));  % int y^p/(y-x)
    for ii = p:-1:2
        P( ii,ifr) = (P(ii+1,ifr)-c(ii))./x(ifr).'; % backward recursion
    end
end

Q = zeros(p,N); % compute q_k from p_k via Helsing 2009 eqn (18)... (p even!)
% Note a rot ang appears here too...  4/17/18
%gam = exp(1i*pi/4); % 1i;  % moves a branch arc as in p_1
%if side == 'e', gam = conj(gam); end   % note gam is a phase, rots branch cut
Q(1:2:end,:) = P(2:2:end,:) - repmat(log((1-x.').*(-1-x.')),[ceil(p/2) 1]); % guessed!
% (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
Q(2:2:end,:) = P(3:2:end,:) - repmat(log(gam) + log((1-x.')./(gam*(-1-x.'))),[floor(p/2) 1]);  % same cut as for p_1
% Seems like abs fails - we must be using complex SLP ? :
%Q(1:2:end,:) = P(2:2:end,:) - repmat(log(abs(1-x.'))+log(abs(-1-x.')),[p/2 1]);
% (-1)^k, k odd, note each log has branch cut in semicircle from -1 to 1
%Q(2:2:end,:) = P(3:2:end,:) - repmat(log(abs(1-x.'))-log(abs(-1-x.')),[p/2 1]);
Q = Q.*repmat(1./(1:p)',[1 N]); % k even
warning('off','MATLAB:nearlySingularMatrix'); % solve for special weights...
A = real((V.'\Q).'.*repmat((1i*s.nx)',[N 1])*zsc)/(2*pi*abs(zsc));
A = A*abs(zsc) - log(abs(zsc))/(2*pi)*repmat(abs(s.wxp)',[N 1]); % unscale, yuk
if nargout>1
    P = P(1:end-1,:);  % trim P back to p rows since kernel is like DLP
    Az = (V.'\P).'*(1/(2*pi)).*repmat((1i*s.nx)',[N 1]); % solve spec wei
    if nargout == 2
        A1 = Az;
    else
        A1 = real(Az); A2 = -imag(Az);
    end
end
end

function P = interpmat(n,m, qntype) % interpolation matrix from n-pt to m-pt Gauss nodes
% INTERPMAT - create interpolation matrix from n-pt to m-pt Gauss nodes
%
% P = interpmat(n,m) returns a m*n matrix which maps func values on n-pt Gauss-
% Legendre nodes on [-1,1] to values on m-pt nodes.
% Does it the Helsing way via backwards-stable ill-cond Vandermonde solve.
if m==n, P = eye(n); return, end
if qntype=='G', x = gauss(n); y = gauss(m); 
else, x = cheby(n); y = cheby(m); end 
V = ones(n); for j=2:n, V(:,j) = V(:,j-1).*x; end % Vandermonde, original nodes
R = ones(m,n); for j=2:n, R(:,j) = R(:,j-1).*y; end % monomial eval matrix @ y
P = (V'\R')';                                       % backwards-stable solve
end

