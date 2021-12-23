function tIC = chebtIC(xgrid,ygrid,xdlim,ydlim,s)
% chebfun to find t intersection....
% 
% Hai 05/08/21

addpath('/home/hszhu/Documents/chebfun');
t = chebfun('t',[0,2*pi]);
if nargin < 5
    a = (3-sqrt(3))/(3+sqrt(3)); w = 5; s = 1/(1+a)*(1 + a*cos(w*t))*exp(1i*t);
end
if ischar(s)    % take string input from python...
    eval(append(s,';'))
end

tol = 1e-10;
t_x = []; 
for k=1:numel(xgrid)
    ylo = ydlim(1); yhi = ydlim(2);  % control verical lines' limits  
    s_xgrid = xgrid(k) + 1i*chebfun('t',[ylo,yhi]); % chebfun version of vertical lines at xgrid(k)
    
    dx = [ylo yhi 0 2*pi];    % specify chebfun domain
    fx = chebfun2(@(t2,t1) s(t1)-s_xgrid(t2),dx);  % rootfinding
    rx = roots(real(fx),imag(fx));              % calculate intersections
    if ~isempty(rx)
        t_x = [t_x;rx(:,2)];
    end
end

t_y = [];
for k=1:numel(ygrid)
    xlo = xdlim(1); xhi = xdlim(2); 
    s_ygrid = chebfun('t',[xlo,xhi]) + 1i*ygrid(k);
    
    dx = [xlo xhi 0 2*pi];    % specify chebfun domain
    fy = chebfun2(@(t2,t1) s(t1)-s_ygrid(t2),dx);  % rootfinding
    ry = roots(real(fy),imag(fy));              % calculate intersections
    if ~isempty(ry)
        t_y = [t_y;ry(:,2)];
    end
end
tIC = uniquetol(mod([t_x;t_y],2*pi),100*tol); % combined interceptions of s with regular mesh, within some tolerance


end