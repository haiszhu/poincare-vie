function ctricell = cell_strcut(s,BaseT,EndX)
% only a template for more general geometry....
% 
% curved triangle 
%
% add template to support boundary nodes
%
% Hai 12/22/21, 

ctricell = [];
for k=1:4
    % curved triangle
    quadBaseT1 = BaseT(k); quadBaseT2 = BaseT(k+1);
    x1 = real(s.Z(quadBaseT1)); x2 = real(s.Z(quadBaseT2)); x3 = real(EndX(k));
    y1 = imag(s.Z(quadBaseT1)); y2 = imag(s.Z(quadBaseT2)); y3 = imag(EndX(k));
    lambda   = @(xi) real(s.Z(quadBaseT1 - xi*(quadBaseT1-quadBaseT2)));
    mu       = @(xi) imag(s.Z(quadBaseT1 - xi*(quadBaseT1-quadBaseT2)));
    lambda_d = @(xi) (quadBaseT2-quadBaseT1)*real(s.Zp(quadBaseT1 - xi*(quadBaseT1-quadBaseT2)));
    mu_d     = @(xi) (quadBaseT2-quadBaseT1)*imag(s.Zp(quadBaseT1 - xi*(quadBaseT1-quadBaseT2)));
    x        = @(xi,eta) (1-xi-eta)*x1 + xi*x2 + eta*x3 + (1-xi-eta)./(1-xi).*(lambda(xi)-(1-xi)*x1-xi*x2);
    y        = @(xi,eta) (1-xi-eta)*y1 + xi*y2 + eta*y3 + (1-xi-eta)./(1-xi).*(    mu(xi)-(1-xi)*y1-xi*y2);
    dxdxi    = @(xi,eta) -x1 + x2 + (1-xi-eta)./(1-xi).*(lambda_d(xi)+x1-x2) + (lambda(xi)-(1-xi)*x1-xi*x2).*(-1./(1-xi)+(1-xi-eta)./(1-xi).^2);
    dxdeta   = @(xi,eta) -x1 + x3 - (lambda(xi)-(1-xi)*x1-xi*x2)./(1-xi);
    dydxi    = @(xi,eta) -y1 + y2 + (1-xi-eta)./(1-xi).*(    mu_d(xi)+y1-y2) + (    mu(xi)-(1-xi)*y1-xi*y2).*(-1./(1-xi)+(1-xi-eta)./(1-xi).^2);
    dydeta   = @(xi,eta) -y1 + y3 - (    mu(xi)-(1-xi)*y1-xi*y2)./(1-xi);
    ctricell{k}.Jac = @(xi,eta) [dxdxi(xi, eta), dxdeta(xi, eta); dydxi(xi, eta), dydeta(xi, eta)];
    ctricell{k}.Jacdet = @(xi,eta) abs(dxdxi(xi, eta).*dydeta(xi, eta) - dxdeta(xi, eta).*dydxi(xi, eta));
    ctricell{k}.xmap = x; ctricell{k}.ymap = y;
    ctricell{k}.vert = [x1,x2,x3;y1,y2,y3];
    if isfield(s,'t')
        idx = (s.t>=quadBaseT1)&(s.t<=quadBaseT2);
        ctricell{k}.bdidx = idx;
        ctricell{k}.bdxi = (s.t(idx)-quadBaseT1)/(quadBaseT2-quadBaseT1);
        ctricell{k}.Zt = 2*pi/3*(s.t(idx)-quadBaseT1)/(quadBaseT2-quadBaseT1) + 2*pi/3;
    end

    % boundary
    triNodes = [x3,x1,x2] + 1i * [y3,y1,y2];
    ctricell{k}.Z  = @(t) ((t>=0)&(t<2*pi/3)).*(triNodes(1)+t/(2*pi/3)*(triNodes(2)-triNodes(1)))...
                        + ((t>=2*pi/3)&(t<4*pi/3)).*s.Z(quadBaseT1 - (t-2*pi/3)/(2*pi/3)*(quadBaseT1-quadBaseT2))...
                        + ((t>=4*pi/3)&(t<=2*pi)).*(triNodes(3)+(t-4*pi/3)/(2*pi/3)*(triNodes(1)-triNodes(3)));    

    % form straight triangle
    lambda = @(xi) x1 + xi*(x2-x1);
    mu = @(xi) y1 + xi*(y2-y1);
    x = @(xi,eta) (1-xi-eta)*x1 + xi*x2 + eta*x3 + (1-xi-eta)./(1-xi).*(lambda(xi)-(1-xi)*x1-xi*x2);
    y = @(xi,eta) (1-xi-eta)*y1 + xi*y2 + eta*y3 + (1-xi-eta)./(1-xi).*(mu(xi)-(1-xi)*y1-xi*y2);
    ctricell{k}.xmap0 = x; ctricell{k}.ymap0 = y;
    
    % boundary
    triNodes0 = [0,0,1] + 1i * [1,0,0];
    ctricell{k}.Z0 = @(t) ((t>=0)&(t<2*pi/3)).*(triNodes0(1)+t/(2*pi/3)*(triNodes0(2)-triNodes0(1)))...
                        + ((t>=2*pi/3)&(t<4*pi/3)).*(triNodes0(2)+(t-2*pi/3)/(2*pi/3)*(triNodes0(3)-triNodes0(2)))...
                        + ((t>=4*pi/3)&(t<=2*pi)).*(triNodes0(3)+(t-4*pi/3)/(2*pi/3)*(triNodes0(1)-triNodes0(3)));
    
%     keyboard
end
end

