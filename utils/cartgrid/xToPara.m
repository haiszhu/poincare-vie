function xToP = xToPara(s)
% from position x to parametric, for an irregular cell
%
% s needs to have order p, nodes(so num of sides of this Ir cell is known),
% and apparently s.Z
%
% tolerance is set to 1e-13...
%
%
% Hai 06/02/21

if nargin == 0, test_xToP; return, end

ref = 3; % default num of panels for each side
p = ceil(3/2*s.p); % use a higher order to find cusp

% compute centroid anyway (might be unnecessary)
nside = numel(s.node);  
s.p = p; s.tpan = linspace(0.0,2*pi,ref*nside+1); 
[s, ~, ~] = quadr(s, [], 'p', 'G');
sIrwst = s.ws.*s.tang; sIrarea = (-1/2*(real(sIrwst))'*(imag(s.x)) + 1/2*(imag(sIrwst))'*(real(s.x)));
sIrCenx = (-1/3*(real(sIrwst))'*(real(s.x).*imag(s.x)) + 1/3*(imag(sIrwst))'*(real(s.x).^2))/sIrarea; sIrCeny = (-1/3*(real(sIrwst))'*(imag(s.x).^2) + 1/3*(imag(sIrwst))'*(imag(s.x).*real(s.x)))/sIrarea;
sIrsc = sIrCenx+1i*sIrCeny; 

[x0,~] = gauss(p);
baseAngle = mod(angle(s.node(1)-sIrsc),2*pi); % the angle of 1st node point...
tmpAngle = unique([0,mod(angle((s.node(2:end)-sIrsc)*exp(-baseAngle*1i)),2*pi),2*pi]); 
tmp = linspace(0,1-1/ref,ref)';
sectionAngle = (1-tmp)*tmpAngle(1:end-1)+tmp*tmpAngle(2:end);
sectionAngle = [sectionAngle(:)',2*pi];
nodeAngle = (1-(x0+1)/2)*sectionAngle(1:end-1) + (x0+1)/2*sectionAngle(2:end);
nodeAngle = nodeAngle(:)';

eps = 1e-11; dt = 0.001;
% starts from 0, with step size dt
t0 = 0+eps; % unless one side is really tiny, this is unlikely to be a bad starting point... (otherwise the mod, angle, function sometimes return 2*pi for t0=0)
T = [0,nodeAngle(1)]; tIC = NaN(size(nodeAngle)); count = 0;
while t0<2*pi && count+1<=numel(nodeAngle)
    while (mod(angle((s.Z(t0)-sIrsc)*exp(-baseAngle*1i)),2*pi)>=T(1)) && (mod(angle((s.Z(t0)-sIrsc)*exp(-baseAngle*1i)),2*pi)<T(2)) % move on until get into next node...
        t0 = t0 + dt;
    end
    t0 = t0 - dt; tl = t0; tr = t0+dt; % there is an intersection within [tl,tr], let't bisection... 
    tmid = 1/2*(tl+tr);
    % try to apply a complex rotation here.
    while abs(mod(angle((s.Z(tmid)-sIrsc)*exp(-baseAngle*1i)),2*pi)-T(2))>eps && abs(tr-tl)>eps
        if mod(angle((s.Z(tmid)-sIrsc)*exp(-baseAngle*1i)),2*pi)-T(2)>0 
            tr = tmid;
        else
            tl = tmid;
        end
        tmid = 1/2*(tl+tr);
    end
    count = count+1;
    tIC(count) = tmid;
    t0 = tIC(count)+dt;
    if count+1<=numel(nodeAngle)
        T = nodeAngle(count:count+1);
    end
end

% find panel end
t0 = 0+eps; 
T = [sectionAngle(1),sectionAngle(2)]; tICplo = NaN(size(sectionAngle)); tICplo(1) = 0; count = 1;
while t0<2*pi && count+1<numel(sectionAngle)
    while (mod(angle((s.Z(t0)-sIrsc)*exp(-baseAngle*1i)),2*pi)>=T(1)) && (mod(angle((s.Z(t0)-sIrsc)*exp(-baseAngle*1i)),2*pi)<T(2)) % move on until get into next node...
        t0 = t0 + dt;
    end
    t0 = t0 - dt; tl = t0; tr = t0+dt; % there is an intersection within [tl,tr], let't bisection... 
    tmid = 1/2*(tl+tr);
    % try to apply a complex rotation here.
    while abs(mod(angle((s.Z(tmid)-sIrsc)*exp(-baseAngle*1i)),2*pi)-T(2))>eps && abs(tr-tl)>eps
        if mod(angle((s.Z(tmid)-sIrsc)*exp(-baseAngle*1i)),2*pi)-T(2)>0 
            tr = tmid;
        else
            tl = tmid;
        end
        tmid = 1/2*(tl+tr);
    end
    count = count+1;
    tICplo(count) = tmid;
    t0 = tICplo(count)+dt;
    T = sectionAngle(count:count+1);
end

% build phase to parametric map...
tmp = mod(angle((s.Z(tICplo)-sIrsc)*exp(-baseAngle*1i)),2*pi); tmp(1) = 0; tmp(end) = 2*pi;
V = ones(p,p); for j=2:p, V(:,j) = V(:,j-1).*x0; end 
xToP =@(x) xToParaFun(x,p,V,tIC,tmp,sIrsc,baseAngle);

end

function t = xToParaFun(x0,p,V,tIC,tmp,sIrsc,baseAngle)

tol = 1e-13;
x = x0(:); tmp = tmp(:)';
xAngle = mod(angle((x-sIrsc)*exp(-baseAngle*1i)),2*pi); % angles for column vec x
idx = sum(xAngle>tmp,2); idx(idx==0) = 1; % figure out which panel
xAnglemid = (tmp(idx)+tmp(idx+1))/2; xAnglemid = xAnglemid(:);  % column vec of panel mid point
xAngleplen = tmp(idx+1)-tmp(idx); xAngleplen = xAngleplen(:);   % column vec of panel length
relAngle = 2*(xAngle-xAnglemid)./xAngleplen;    % standard panel for interpolation
R = ones(numel(x),p); for j=2:p, R(:,j) = R(:,j-1).*relAngle; end % target interpolation matrix (V is source)
t = sum((V'\R')'.*tIC((idx-1)*p+(1:p)),2); % tIC((idx-1)*p+(1:p)) is a numel(x) X p matrix, parametric info for interpolation
                                           % (V'\R')' has each row as the interpolation matrix for x(i)
t = reshape(t,size(x0)); t(abs(t-2*pi)<tol) = 0;
end

function test_xToP()

p = 10;

s.Z = @(t)  sqrt(2)*(  ((t>=0)&(t<pi/2)).*(1./(sin(t)+cos(t))).*(cos(t)+1i*sin(t))*exp(1i*pi/4)...
                     + ((t>=pi/2)&(t<pi)).*(1./(sin(t-pi/2)+cos(t-pi/2))).*(cos(t-pi/2)+1i*sin(t-pi/2))*exp(1i*3*pi/4)...
                     + ((t>=pi)&(t<3*pi/2)).*(1./(sin(t-pi)+cos(t-pi))).*(cos(t-pi)+1i*sin(t-pi))*exp(1i*5*pi/4)...
                     + ((t>=3*pi/2)&(t<=2*pi)).*(1./(sin(t-3*pi/2)+cos(t-3*pi/2))).*(cos(t-3*pi/2)+1i*sin(t-3*pi/2))*exp(1i*7*pi/4));
s.node = s.Z(linspace(0,2*pi,5)); s.node = s.node(1:4);
s.sc = 0; s.p = p;
xToP = xToPara(s);

[XX,YY] = meshgrid(cheby(12));
t.x = XX+1i*YY;
test1 = xToP(t.x);
test2 = mod(mod(angle(t.x-s.sc)-pi/4,2*pi),2*pi);

diff = abs(test1-test2);

keyboard

end