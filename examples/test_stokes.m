% test Stokes volume potential eval on a star shape geometry
%
%
% Hai 12/22/21

addpath('../')
setup()

% geometry
a = -.1; w = 4;           % smooth wobbly radial shape params...
R = @(t) (1 + a*cos(w*t))*1; Rp = @(t) -w*a*sin(w*t); Rpp = @(t) -w*w*a*cos(w*t);
s.Z = @(t) R(t).*exp(1i*t); s.Zp = @(t) (Rp(t) + 1i*R(t)).*exp(1i*t);
s.Zpp = @(t) (Rpp(t) + 2i*Rp(t) - R(t)).*exp(1i*t);

% set up four cells (each quadrant is one irregular cell)
BaseT = linspace(0,2*pi,5); % curved triangular patches are based off s.Z, counterclockwise ?, how about exterior
EndX = (0+1i*0)*ones(1,numel(BaseT)-1); % the other node that has two straight sides
s.tpan = linspace(0,2*pi,21); s.p = 10; s = quadr(s,[],'p','G'); s_out = s;
ctricell = cell_strcut(s,BaseT,EndX); 

% boundary targets...
for k=1:numel(ctricell)
    sk = ctricell{k};
    sxbdk = sk.xmap(sk.bdxi,0*sk.bdxi) + 1i*sk.ymap(sk.bdxi,0*sk.bdxi);
    Ztbd = sk.Zt; % tmp = sk.Z(Ztbd);
end

% for reference triangle (which is the same for every element)
vert0 = [[1,0,0];[0,1,0]]';
x1 = vert0(1,1); x2 = vert0(2,1); x3 = vert0(3,1);
y1 = vert0(1,2); y2 = vert0(2,2); y3 = vert0(3,2);
lambda = @(xi) x1 + xi*(x2-x1);
mu = @(xi) y1 + xi*(y2-y1);
xmap = @(xi,eta) (1-xi-eta)*x1 + xi*x2 + eta*x3 + (1-xi-eta)./(1-xi).*(lambda(xi)-(1-xi)*x1-xi*x2);
ymap = @(xi,eta) (1-xi-eta)*y1 + xi*y2 + eta*y3 + (1-xi-eta)./(1-xi).*(mu(xi)-(1-xi)*y1-xi*y2);

% target
norder = 15;
[uvs,wts]=get_vioreanu_nodes(norder);
T.x = [];

% volume density
mu1 = @(x,y) 2*y.*(1+cos(x.*y)); mu2 = @(x,y) -2*x.*(1-cos(x.*y));

% alpert
p = 16;
[atl, awl] = QuadNodesInterval(0, 1, 16, 0, 3, 1, p); % alpert with left end correction
t0a = atl(:); w0a = awl(:);

% 
val1 = []; val2 = [];
for k=1:4 % how to deal with close targets? 
    tbd.x = s.x(((k-1)*end/4+1):k*end/4);
    t.x = ctricell{k}.xmap(uvs(1,:),uvs(2,:)) + 1i*ctricell{k}.ymap(uvs(1,:),uvs(2,:)); t.x = t.x(:);
    
    % volume potential
    val1{k} = zeros(size([t.x(:);tbd.x(:)])); val2{k} = val1{k};
    % self patch
    sk = ctricell{k}; % self patch
    Zt = findt(uvs(1,:)+1i*uvs(2,:)); Ztbd = ctricell{k}.bdxi;
    % combine interior and boundary targets
    R0 = [uvs,[ctricell{k}.bdxi(:),0*ctricell{k}.bdxi(:)]'];
    t.x = [t.x;tbd.x]; T.x = [T.x;t.x];
    Zt = [Zt(:);Ztbd(:)];

    % loop over targets
    sk.Z = sk.Z0; % change back to reference space
    for j=1:numel(t.x)
        
        % build target dependent singular quadr
    %     r0 = [uvs(1,j);uvs(2,j)]; % for self, singularity is exact
        r0 = R0(:,j);
        t0 = Zt(j);
        sk.p = 3/2*p; tpan0 = linspace(0,2*pi,10); tlo = tpan0(sum(t0>=tpan0)); thi = tpan0(sum(t0>=tpan0)+1); 
        plen = thi - tlo;
        tpantmp = unique(mod([tpan0,t0-plen*(1/4).^(1:6),t0+plen*(1/4).^(1:6),t0],2*pi));
        tpanIdx = abs([tpantmp(2:end),2*pi]-tpantmp)>1e-13;
        sk.tpan = [tpantmp(tpanIdx),2*pi];
        sk = quadr(sk,[],'p','G');
        r = [real(sk.x),imag(sk.x)]';
        % in reference space
        tr_x0 = t0a(:)*(r(1,:)-r0(1))+r0(1); tr_y0 = t0a(:)*(r(2,:)-r0(2))+r0(2); tw = w0a(:).*t0a(:); % figure(),plot(tr_x,tr_y,'.'); axis equal
        Dx = sk.ws.*real(sk.tang); Dy = sk.ws.*imag(sk.tang);
        tr_ws0 = tw*(-(r(2,:)-r0(2)).*Dx(:)'+(r(1,:)-r0(1)).*Dy(:)'); 
        % in physical space
        Jacdet = sk.Jacdet(tr_x0(:),tr_y0(:));
        tr_x = sk.xmap(tr_x0(:),tr_y0(:)); tr_y = sk.ymap(tr_x0(:),tr_y0(:));
        tr_ws = tr_ws0(:).*Jacdet(:); 
        
        % kernel computation
        x = tr_x(:)'; y = tr_y(:)'; wtmp1 = tr_ws(:).*mu1(x(:),y(:)); wtmp2 = tr_ws(:).*mu2(x(:),y(:)); 
        vecm = (real(t.x(j))-x)+1i*(imag(t.x(j))-y);
        logPart = -1/(4*pi)*log(abs(vecm));
        crsPart =  1/(4*pi)*(real(vecm).*imag(vecm)./abs(vecm).^2);
        StoSLVPmatm11 = (logPart + 1/(4*pi)*(real(vecm).^2./abs(vecm).^2))*wtmp1;
        StoSLVPmatm12 = crsPart*wtmp2; StoSLVPmatm21 = crsPart*wtmp1;
        StoSLVPmatm22 = (logPart + 1/(4*pi)*(imag(vecm).^2./abs(vecm).^2))*wtmp2;
        val1{k}(j) = val1{k}(j) + StoSLVPmatm11 + StoSLVPmatm12;
        val2{k}(j) = val2{k}(j) + StoSLVPmatm21 + StoSLVPmatm22;
    end

    % neighbor patches
    for l=[1:k-1,k+1:4] % things get tricky, maybe assume curved side effect is small
    
        sk = ctricell{l}; % loop over patches
    
        % build inverse map roughly
        vert = sk.vert;
        x1 = vert(1,1); x2 = vert(1,2); x3 = vert(1,3);
        y1 = vert(2,1); y2 = vert(2,2); y3 = vert(2,3);
        xietaMap = 1/((x2-x1)*(y3-y1)-(y2-y1)*(x3-x1))*[y3-y1,x1-x3;y1-y2,x2-x1];
        % now for the real targets... that might be close
        xieta = xietaMap*[real(t.x(:))-x1,imag(t.x(:))-y1]';
    
        % loop over targets
        sk.Z = sk.Z0; % change back to reference space
        for j=1:numel(t.x)
    
            % build target dependent near singular quadr
            targj = [real(t.x(j));imag(t.x(j))]; % target in physical space
            % 
            xieta0 = xieta(:,j); % (xi,eta) outside reference triangle
            r0 = [xmap(xieta0(1),xieta0(2)); ymap(xieta0(1),xieta0(2))]; % same map, becasue findt is dumb
            t0 = findt(xieta0(1)+1i*xieta0(2)); % 
            rstar = ctricell{k}.Z(t0); 
            % build panel boundary
            sk.p = 3/2*p; tpan0 = linspace(0,2*pi,10); tlo = tpan0(sum(t0>=tpan0)); thi = tpan0(sum(t0>=tpan0)+1); 
            plen = thi - tlo;
            tpantmp = unique(mod([tpan0,t0-plen*(1/4).^(1:6),t0+plen*(1/4).^(1:6),t0],2*pi));
            tpanIdx = abs([tpantmp(2:end),2*pi]-tpantmp)>1e-13;
            sk.tpan = [tpantmp(tpanIdx),2*pi];
            sk = quadr(sk,[],'p','G');
            r = [real(sk.x),imag(sk.x)]';
            % k-r in reference space
            XiEta1 = sk.Z(t0); xieta1 = [real(XiEta1);imag(XiEta1)];
            a = 0; b = 1; dist = -abs((norm(xieta1-[1/3;1/3])-norm(xieta0-[1/3;1/3]))/norm(xieta1-[1/3;1/3]));
            [xx, ww] = neighQuad(a, b, dist, 'left'); krtl = xx/(b-a); krwl = ww/(b-a);
            t0kr = krtl(:); w0kr = krwl(:);
            tr_x0 = t0kr(:)*(r(1,:)-xieta1(1))+xieta1(1); tr_y0 = t0kr(:)*(r(2,:)-xieta1(2))+xieta1(2); tw = w0kr(:).*t0kr(:);
           
            Dx = sk.ws.*real(sk.tang); Dy = sk.ws.*imag(sk.tang);
            tr_ws0 = tw*(-(r(2,:)-xieta1(2)).*Dx(:)'+(r(1,:)-xieta1(1)).*Dy(:)'); 
            % in physical space
            Jacdet = sk.Jacdet(tr_x0(:),tr_y0(:));
            tr_x = sk.xmap(tr_x0(:),tr_y0(:)); tr_y = sk.ymap(tr_x0(:),tr_y0(:));
            tr_ws = tr_ws0(:).*Jacdet(:); 
    
            % kernel computation
            x = tr_x(:)'; y = tr_y(:)'; wtmp1 = tr_ws(:).*mu1(x(:),y(:)); wtmp2 = tr_ws(:).*mu2(x(:),y(:)); 
            vecm = (real(t.x(j))-x)+1i*(imag(t.x(j))-y);
            logPart = -1/(4*pi)*log(abs(vecm));
            crsPart =  1/(4*pi)*(real(vecm).*imag(vecm)./abs(vecm).^2);
            StoSLVPmatm11 = (logPart + 1/(4*pi)*(real(vecm).^2./abs(vecm).^2))*wtmp1;
            StoSLVPmatm12 = crsPart*wtmp2; StoSLVPmatm21 = crsPart*wtmp1;
            StoSLVPmatm22 = (logPart + 1/(4*pi)*(imag(vecm).^2./abs(vecm).^2))*wtmp2;
            val1{k}(j) = val1{k}(j) + StoSLVPmatm11 + StoSLVPmatm12;
            val2{k}(j) = val2{k}(j) + StoSLVPmatm21 + StoSLVPmatm22;
           
    %         keyboard
            
        end
    
    end

end

%% a global approach
% old way, direct computation without forming the basis
t.x = T.x(:);
val01 = zeros(size(t.x)); val02 = val01;
tref = linspace(0, 2*pi, 4001); tref = tref(1:end-1); sxtmp = s.Z(tref);
[~,t0idxc] = min(abs(t.x(:)-sxtmp),[],2);
Zt2 = tref(t0idxc);
for k=1:numel(t.x)
    r0 = [real(t.x(k));imag(t.x(k))];
    t0 = Zt2(k);
    s.p = 3/2*p; tpan0 = linspace(0,2*pi,33); tlo = tpan0(sum(t0>=tpan0)); thi = tpan0(sum(t0>=tpan0)+1); 
    plen = thi - tlo;
    tpantmp = unique(mod([tpan0,t0-plen*(1/4).^(1:6),t0+plen*(1/4).^(1:6),t0],2*pi));
    tpanIdx = abs([tpantmp(2:end),2*pi]-tpantmp)>1e-13;
    s.tpan = [tpantmp(tpanIdx),2*pi];
    s = quadr(s,[],'p','G');
    r = [real(s.x),imag(s.x)]';

    %
    tr_x = t0a(:)*(r(1,:)-r0(1))+r0(1); tr_y = t0a(:)*(r(2,:)-r0(2))+r0(2); tw = w0a(:).*t0a(:); % figure(),plot(tr_x,tr_y,'.'); axis equal
    Dx = s.ws.*real(s.tang); Dy = s.ws.*imag(s.tang);
    tr_ws = tw*(-(r(2,:)-r0(2)).*Dx(:)'+(r(1,:)-r0(1)).*Dy(:)'); 

    % volume integral of Stokes
    x = tr_x(:)'; y = tr_y(:)'; wtmp1 = tr_ws(:).*mu1(x(:),y(:)); wtmp2 = tr_ws(:).*mu2(x(:),y(:)); 
    vecm = (r0(1)-x)+1i*(r0(2)-y);
    logPart = -1/(4*pi)*log(abs(vecm));
    crsPart =  1/(4*pi)*(real(vecm).*imag(vecm)./abs(vecm).^2);
    StoSLVPmatm11 = (logPart + 1/(4*pi)*(real(vecm).^2./abs(vecm).^2))*wtmp1;
    StoSLVPmatm12 = crsPart*wtmp2; StoSLVPmatm21 = crsPart*wtmp1;
    StoSLVPmatm22 = (logPart + 1/(4*pi)*(imag(vecm).^2./abs(vecm).^2))*wtmp2;

    val01(k) = StoSLVPmatm11 + StoSLVPmatm12;
    val02(k) = StoSLVPmatm21 + StoSLVPmatm22;

end

% volume potential error (estimated)
val1 = vertcat(val1{:}); val2 = vertcat(val2{:});
err1 = reshape(abs(val1-val01),[],4); err2 = reshape(abs(val2-val02),[],4);

% volume potential computed in reference space
Val1 = reshape(val1,[],4); Val2 = reshape(val2,[],4);

% volume integral done... solve Stokes problem
ubd = -real(s_out.x).^2.*imag(s_out.x) + 1i*real(s_out.x).*imag(s_out.x).^2;
vValbd = Val1((end-numel(s_out.x)/4+1):end,:) + 1i*Val2((end-numel(s_out.x)/4+1):end,:);
vValbd = vValbd(:);
uHombd = ubd - vValbd;
Asto = StoSLPSpecialMat(s_out,s_out,'i');
tau = Asto\[real(uHombd);imag(uHombd)];

% compute homogeneous soln
t.x = reshape(T.x,[],4); t.x = t.x(1:(end-numel(s_out.x)/4),:); t.x = t.x(:);
uHom = StoSLPSpecialMat(t,s_out,'i')*tau; uHom = uHom(1:end/2) + 1i*uHom(end/2+1:end);
uVal = Val1(1:(end-numel(s_out.x)/4),:) + 1i*Val2(1:(end-numel(s_out.x)/4),:); uVal = uVal(:);
uNum = uHom + uVal;
uExa = -real(t.x).^2.*imag(t.x) + 1i*real(t.x).*imag(t.x).^2;

% error
Err = abs(uNum-uExa);

figure(1),clf,scatter(real(t.x),imag(t.x),[],log10(Err)), title('forced Stokes error'),colorbar

keyboard

