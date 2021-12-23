function [ftval,cMat] = OrthoApprox(sx,ws,p,fval,tx,cMat)
% Approx with orthogonal polynomials on irregular domain...
%
% quadrature nodes -- sx, as long as num(sx) within a few hundreds, this should be reasonably fast
% quadrature weights -- ws
% approximation order -- p
% sampled function value -- fval at sx
% target points -- tx
% Gram-Schmidt orthogonalization process - cMat (don't have to provide this)
%
% Hai 10/13/21

if nargin==0, test_OrthoApprox; return; end

% smooth quadr nodes and weights
x = real(sx(:)); y = imag(sx(:)); w = ws(:);

% start orthogonalization process

% at quadr
lMat0 = NaN(numel(x),p^2); % this records legendre polynomial evaluated at quadrature nodes..
lMat0(:,1) = [ones(size(x))]; % start with constant 
if nargin < 6 % without cMat, we have to find orthogonal polynomials...
    cMat = eye(p^2); % process along the way, this records the projection coefficients;
                 % the 1st few also need to be changed for irregular boxes
    % initialize new basis (multiple by the corresponding x or y with proper coefficient)
    for j = 1:p-1
        % at quadr
        for k=1:j-1
            lMat0(:,j^2+2*k-1) = bsxfun(@times,(2*j-1)/j*x,lMat0(:,(j-1)^2+2*k-1));
            lMat0(:,j^2+2*k) = bsxfun(@times,(2*j-1)/j*y,lMat0(:,(j-1)^2+2*k));
        end
        lMat0(:,j^2+2*j-1) = bsxfun(@times,(2*j-1)/j*x,lMat0(:,(j-1)^2+2*j-1));
        lMat0(:,j^2+2*j) = bsxfun(@times,(2*j-1)/j*y,lMat0(:,(j-1)^2+2*j-1));
        lMat0(:,j^2+2*j+1) = bsxfun(@times,((2*j-1)/j)^2*x.*y,lMat0(:,(j-1)^2+2*j-1));

        % discrete orthogonalization process
        for l = 1:2*j+1
            for k=1:j^2+l-1
                cMat(k,j^2+l) = -sum(lMat0(:,j^2+l).*lMat0(:,k).*w)/sum(lMat0(:,k).^2.*w);
            end
            lMat0(:,j^2+l) = lMat0(:,1:j^2+l)*cMat(1:j^2+l,j^2+l);
        end
    end
else
    % initialize new basis (multiple by the corresponding x or y with proper coefficient)
    for j = 1:p-1
        % at quadr
        for k=1:j-1
            lMat0(:,j^2+2*k-1) = bsxfun(@times,(2*j-1)/j*x,lMat0(:,(j-1)^2+2*k-1));
            lMat0(:,j^2+2*k) = bsxfun(@times,(2*j-1)/j*y,lMat0(:,(j-1)^2+2*k));
        end
        lMat0(:,j^2+2*j-1) = bsxfun(@times,(2*j-1)/j*x,lMat0(:,(j-1)^2+2*j-1));
        lMat0(:,j^2+2*j) = bsxfun(@times,(2*j-1)/j*y,lMat0(:,(j-1)^2+2*j-1));
        lMat0(:,j^2+2*j+1) = bsxfun(@times,((2*j-1)/j)^2*x.*y,lMat0(:,(j-1)^2+2*j-1));

        % discrete orthogonalization process
        for l = 1:2*j+1
            lMat0(:,j^2+l) = lMat0(:,1:j^2+l)*cMat(1:j^2+l,j^2+l);
        end
    end
end

% now with orthogonal polynomials at these quadr nodes 
approxMat = zeros(p^2,1); rhs = zeros(p^2,1); fval = fval(:);
for j=1:p^2
    approxMat(j) = sum(lMat0(:,j).^2.*w);
    rhs(j) = sum(lMat0(:,j).*fval.*w);
end
coeff = rhs./approxMat;

% targets
targx = real(tx(:)); targy = imag(tx(:));
lMat = NaN(numel(targx),p^2); 
lMat(:,1) = [ones(size(targx))];
for j = 1:p-1
    % recursively moving forwards to higher order...
    for k=1:j-1
        lMat(:,j^2+2*k-1) = bsxfun(@times,(2*j-1)/j*targx,lMat(:,(j-1)^2+2*k-1));
        lMat(:,j^2+2*k) = bsxfun(@times,(2*j-1)/j*targy,lMat(:,(j-1)^2+2*k));
    end
    lMat(:,j^2+2*j-1) = bsxfun(@times,(2*j-1)/j*targx,lMat(:,(j-1)^2+2*j-1));
    lMat(:,j^2+2*j) = bsxfun(@times,(2*j-1)/j*targy,lMat(:,(j-1)^2+2*j-1));
    lMat(:,j^2+2*j+1) = bsxfun(@times,((2*j-1)/j)^2*targx.*targy,lMat(:,(j-1)^2+2*j-1));
    % discrete orthogonalization process
    for l = 1:2*j+1
        lMat(:,j^2+l)  = lMat(:,1:j^2+l)*cMat(1:j^2+l,j^2+l); 
    end
end
ftval = lMat*coeff;

end

function test_OrthoApprox
% simple test ...
%
addpath('../../')
setup()

p = 16; ref = 1; %
f = @(x,y) sin(x) + cos(y) + exp(x.*y);% + sin(3*x) + cos(4*x).*sin(5*y);

if 0
    % smooth volume quadr (need to do more test to see if upsampling is necessary)
    box.Z = @(t) exp(1i*t);
    box.p = ref*p; box.tpan = linspace(0,2*pi,9);                 
    [box, ~, ~] = quadr(box, [], 'p', 'G'); % this is just a shift of sIr.x, one extra quadr is ok...
    [x0b,w0b] = gauss(ref*p); 
    r = [real(box.x),imag(box.x)]';
    tr_x = (x0b(:)+1)/2*(r(1,:)); tr_y = (x0b(:)+1)/2*(r(2,:));
    Dx = box.ws.*real(box.tang); Dy = box.ws.*imag(box.tang);
    tr_ws = 1/4*(w0b(:).*(x0b(:)+1))*(-(r(2,:)).*Dx(:)'+(r(1,:)).*Dy(:)');
    x = tr_x(:); y = tr_y(:); w = tr_ws(:);
else
    a = -.05; b = 6;  
    box.Z = @(t) ((t>=0)&(t<2*pi/2)).*(-1-1i + 2*(1 + a*sin(b*(t)/2)).*exp(1i*(t)/2)) + ((t>=2*pi/2)&(t<3*pi/2)).*(-1+1i*(1-4/pi*(t-2*pi/2))) + ((t>=3*pi/2)&(t<=2*pi)).*(-1+4/pi*(t-3*pi/2)-1i); 
    box.p = ref*p; box.tpan = linspace(0,2*pi,5);
    [box, ~, ~] = quadr(box, [], 'p', 'G'); % this is just a shift of sIr.x, one extra quadr is ok...
    [x0b,w0b] = gauss(ref*p); 
    r = [real(box.x),imag(box.x)]';
    tr_x = (x0b(:)+1)/2*(r(1,:)); tr_y = (x0b(:)+1)/2*(r(2,:));
    Dx = box.ws.*real(box.tang); Dy = box.ws.*imag(box.tang);
    tr_ws = 1/4*(w0b(:).*(x0b(:)+1))*(-(r(2,:)).*Dx(:)'+(r(1,:)).*Dy(:)');
    x = tr_x(:); y = tr_y(:); w = tr_ws(:);
end

% targets
tNum = 51;
[tXX,tYY] = meshgrid(linspace(-1,1,tNum));
tx = tXX(:); ty = tYY(:);

% call approximation
[ftval,~] = OrthoApprox(x(:)+1i*y(:),w(:),p,f(x(:),y(:)),tx(:)+1i*ty(:));

% error
err = abs(ftval - f(tx,ty));

figure(1),clf,
errplot = reshape(log10(abs(err)),tNum,tNum);
pcolor(tXX,tYY,errplot);colorbar
hold on,
plot(box.x,'.-r')
title('approximation error (ignore extrapolation part...)')

figure(2),clf,
surf(tXX,tYY,f(tXX,tYY)); axis equal, hold on
plot3(tx,ty,ftval,'.r')
title('surf plot')

end