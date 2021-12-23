% concentric circle example
% 
% 
% Hai 08/18/21

% profile on

clear all
addpath('..'); setup()
addpath('/home/hszhu/Documents/chebfun');

% -\Delta ue(x,y) = mu(x,y)
mu = @(x,y) 6*sin(6*x)+8*cos(8*(y+0.1))+4*(x.^2+y.^2).*sin(4*x.*y)+3*cos(3*x).*sin(3*y); 
ue = @(x,y) 1/6*sin(6*x)+1/8*cos(8*(y+0.1))+1/4*sin(4*x.*y)+1/6*cos(3*x).*sin(3*y);

% volume mesh
if 1
    gridsize = 1/5; p = 10; % p = 10 for error plot, p = 5; if want to plot quadr nodes and weight figures...
    s0.Z = @(t) 1/3*exp(1i*t); se.Z = @(t) 8/9*exp(1i*t); 
    s0_cheb = 1/3*exp(1i*chebfun('t',[0,2*pi])); se_cheb = 8/9*exp(1i*chebfun('t',[0,2*pi])); 
else
    gridsize = 1/10; p = 10;
    a = .1; b = 8; 
    s0.Z = @(t) 5/9/sqrt(2)*(1 + a*cos(b*t)).*exp(1i*t); 
    se.Z = @(t) 10/9/sqrt(2)*(1 + a/2*cos(b*t)).*exp(1i*t); 
    s0_cheb = 5/9/sqrt(2)*(1 + a*cos(b*chebfun('t',[0,2*pi])))*exp(1i*chebfun('t',[0,2*pi])); 
    se_cheb = 10/9/sqrt(2)*(1 + a/2*cos(b*chebfun('t',[0,2*pi])))*exp(1i*chebfun('t',[0,2*pi])); 
    up = 150; 
end
si{1}.Z = s0.Z; 
temp = se.Z(linspace(0,2*pi,500));
sxlo = 1.1*min(real(temp)); sxhi = 1.1*max(real(temp)); sylo = 1.1*min(imag(temp)); syhi = 1.1*max(imag(temp));
xgrid = gridsize*floor(sxlo/gridsize):gridsize:gridsize*ceil(sxhi/gridsize);
ygrid = gridsize*floor(sylo/gridsize):gridsize:gridsize*ceil(syhi/gridsize);
xdlim = [min(xgrid),max(xgrid)]; ydlim = [min(ygrid),max(ygrid)];
gridsize = xgrid(2)-xgrid(1);
tICe = chebtIC(xgrid,ygrid,xdlim,ydlim,se_cheb);
tICi1 = chebtIC(xgrid,ygrid,xdlim,ydlim,s0_cheb);
se.tIC = tICe; si{1}.tIC = tICi1;
[sR, sRIdxMat] = rbox(xgrid,ygrid,se,si);
[sIre, sIrIdxMate] = irbox(xgrid,ygrid,se,se.tIC,'i'); se.tICkeep = true(size(se.tIC)); 
[sIr{1}, sIrIdxMat1] = irbox(xgrid,ygrid,si{1},si{1}.tIC,'e'); si{1}.tICkeep = true(size(si{1}.tIC));
sI = [sIr{1} sIre]; nirreg = numel(sI);
Nbir = numel(sI); Nbr = numel(sR);

% regular box smooth quadr
ref = 1; [xf,wf] = cheby(ref*p); [XXf,YYf] = meshgrid(xf); WWf = (wf*wf');
for j=1:numel(sR)
    sR{j}.tr = [(gridsize/2)*XXf(:)+real(sR{j}.sc),(gridsize/2)*YYf(:)+imag(sR{j}.sc)]'; sR{j}.tws = (gridsize/2)^2*WWf(:)';
end
% irregular box smooth quadr
ref = 1; [xf,wf] = gauss(ref*p); 
for j=1:numel(sI)
    nside = numel(sI{j}.node); sI{j}.sc0 = sI{j}.sc;
    % find centroid, quadr...
    stmp.Z = sI{j}.Z; stmp.tpan = linspace(0.0,2*pi,nside+1); stmp.p = ref*p;
    stmp = quadr(stmp, [], 'p', 'G');
    swst = stmp.ws.*stmp.tang; area = (-1/2*(real(swst))'*(imag(stmp.x)) + 1/2*(imag(swst))'*(real(stmp.x)));
    sCenx = (-1/3*(real(swst))'*(real(stmp.x).*imag(stmp.x)) + 1/3*(imag(swst))'*(real(stmp.x).^2))/area; 
    sCeny = (-1/3*(real(swst))'*(imag(stmp.x).^2) + 1/3*(imag(swst))'*(imag(stmp.x).*real(stmp.x)))/area;
    sIrsc = sCenx+1i*sCeny;
    r = [real(stmp.x)-sCenx,imag(stmp.x)-sCeny]'; 
    tr_x = (xf(:)+1)/2*(r(1,:)); tr_y = (xf(:)+1)/2*(r(2,:));
    Dx = stmp.ws.*real(stmp.tang); Dy = stmp.ws.*imag(stmp.tang);
    tr_ws = 1/4*(wf(:).*(xf(:)+1))*(-(r(2,:)).*Dx(:)'+(r(1,:)).*Dy(:)');
    % additional attribute
    sI{j}.sc = sIrsc;
    sI{j}.tr = [tr_x(:)+sCenx,tr_y(:)+sCeny]'; sI{j}.tws = tr_ws(:)';
end
sV = [sI sR];

figure(1),clf,hold on,axis equal
for k=1:Nbir, plot(sV{k}.Z(linspace(0,2*pi,100)),'-k','LineWidth',1); end
for j=1:Nbr,  plot(sV{j+Nbir}.Z(linspace(0,2*pi,100)),'--b','LineWidth',1); end

% opts.dt = 0.1; memorygraph('start',opts);
% disp('testing memorygraph: please wait 10 secs...')

% up = 200; % 2.6891e-10
up = 250; delta = 0.001;
[uXX,uYY] = meshgrid(linspace(-1+1/up,1-1/up,up));
[IN1,ON1] = inpolygon(uXX(:),uYY(:),(1+delta)*real(s0.Z(linspace(0,2*pi,200))),(1+delta)*imag(s0.Z(linspace(0,2*pi,200))));
[IN2,ON2] = inpolygon(uXX(:),uYY(:),(1-delta)*real(se.Z(linspace(0,2*pi,200))),(1-delta)*imag(se.Z(linspace(0,2*pi,200))));
t.x = (uXX(~IN1&~ON1&IN2&~ON2)+1i*uYY(~IN1&~ON1&IN2&~ON2)); tr0 = [real(t.x),imag(t.x)]';

% inner circle boundary condition
Np1 = 32*2; s1.p = 16; s1.tpan = linspace(0,2*pi,Np1+1)'; s1.Z = @(t) s0.Z(t); s1 = quadr(s1,[],'p','G'); 
Np2 = 32*2; s2.p = 16; s2.tpan = linspace(0,2*pi,Np2+1)'; s2.Z = @(t) se.Z(t); s2 = quadr(s2,[],'p','G'); 

% concatenate targets and source
tbd.x = [s1.x;s2.x];   % targets on the boundary,
tin.x = tr0(1,:)+1i*tr0(2,:); % 
tR = cellfun(@(d)d.tr,sV,'uniformoutput',0); tR = horzcat(tR{:}); 
tWs = cellfun(@(d)d.tws,sV,'uniformoutput',0); tWs = horzcat(tWs{:}); 
muTxy = mu(tR(1,:),tR(2,:));

% naive volume integral
vs.x = tR(1,:)'+1i*tR(2,:)'; vs.ws = tWs(:); % volume source
vt.x = tbd.x(:); ut.x = tin.x.';

% disp('quadr and t'); memorygraph('label','quadr and t');

[~,vValfmm]=LapSLPFMM2D(vs,0,vt,1,muTxy,4);
[~,uValfmm]=LapSLPFMM2D(vs,0,ut,1,muTxy,4);
vVal = vValfmm; uVal = uValfmm;

% disp('fmm'); memorygraph('label','fmm');

% correction from irregular boxes
[QuadI,~] = aLapSLVPQuad(p,sI,[],tbd,tin);
coeff = zeros(p^2*numel(sI),1);
for k=1:numel(sI)
    coeff((k-1)*p^2+(1:p^2)) = QuadI.coeffFun{k}(mu(real(QuadI.sx{k}),imag(QuadI.sx{k})));
end
vVal = vVal + (QuadI.Abd*coeff)'; 
uVal = uVal + (QuadI.Ain*coeff)';

% correction from regular boxes
[~,QuadR] = aLapSLVPQuad(p,[],sR,tbd,tin);
vVal = vVal + (QuadR.Abd*mu(real(QuadR.sx),imag(QuadR.sx)))';
uVal = uVal + (QuadR.Ain*mu(real(QuadR.sx),imag(QuadR.sx)))';

% disp('eval matrix'); memorygraph('label','eval matrix');

%% volume integral done... solve Laplace problem (a naive setup)
fe1 = ue(real(s1.x),imag(s1.x)); 
fe2 = ue(real(s2.x),imag(s2.x)); 

% rhs
rhs = [fe1;fe2]-vVal(:);
% Laplace self interaction matrix
As1 = LapSLPSpecialMat(s1,s1,'e'); As2 = LapSLPSpecialMat(s2,s2,'i'); 
Amat = [As1,LapSLPmat(s1,s2);LapSLPmat(s2,s1),As2];
% solve
tau = Amat\rhs;
% eval at target
uLap = LapSLPSpecialMat(t,s1,'e')*tau(1:numel(s1.x))+LapSLPSpecialMat(t,s2)*tau(numel(s1.x)+1:end);
% poisson solution
uNum = uLap(:)+uVal(:);

% disp('bie part'); memorygraph('label','bie part');

% error 
Err = NaN(size(uXX));
Err(~IN1&~ON1&IN2&~ON2) = abs(uNum-ue(real(t.x)',imag(t.x)')')/max(abs(uNum(:)));

% volume integral value
valtmp = NaN(size(uXX));
valtmp(~IN1&~ON1&IN2&~ON2) = uVal;

figure(2),clf,surf(uXX,uYY,valtmp,'EdgeColor','None'), title('u^P'); hold on, 
plot3(real(s1.x),imag(s1.x),vVal(1:numel(s1.x)),'k','LineWidth',2)
plot3(real(s2.x),imag(s2.x),vVal(numel(s1.x)+1:end),'k','LineWidth',2)
cmp = getPyPlot_cMap('rainbow', [], [], '"/home/hszhu/.pyenv/shims/python"');
colormap(cmp),colorbar

% figure(3),clf,surf(uXX,uYY,Err), title('error')
figure(3),clf,h = pcolor(uXX,uYY,log10(Err)); hold on; title('error')
colorbar; caxis([-14 -7]); set(h, 'EdgeColor', 'none');
plot(s1.x,'k','LineWidth',2); plot(s2.x,'k','LineWidth',2); axis equal
xlim([min(tR(1,:))-1e-2 max(tR(1,:))+1e-2])
ylim([min(tR(2,:))-1e-2 max(tR(2,:))+1e-2])
colormap(cmp),

% memorygraph('label','density solve done');
% 
% [b et ct c las lat] = memorygraph('plot');
% memorygraph('done');

% profile viewer

keyboard

% cmp = getPyPlot_cMap('rainbow', [], [], '"/home/hszhu/.pyenv/shims/python"');
% colormap(cmp),
% 
% figure(1),clf,scatter(tR(1,:),tR(2,:),2,tWs,'filled')
% colorbar, title('quadrature'),axis equal
% cmp = getPyPlot_cMap('rainbow', [], [], '"/home/hszhu/.pyenv/shims/python"');
% colormap(cmp),
% xlim([min(tR(1,:))-1e-2 max(tR(1,:))+1e-2])
% ylim([min(tR(2,:))-1e-2 max(tR(2,:))+1e-2])
% 
% figure(1), print -dpng -r600 ex4_1.png; system('convert -trim ex4_1.png eps2:ex4_1.eps');
% figure(2), print -dpng -r600 ex4_2.png; system('convert -trim ex4_2.png eps2:ex4_2.eps');
% figure(3), print -dpng -r600 ex4_3.png; system('convert -trim ex4_3.png eps2:ex4_3.eps');
% 
% figure(1), print -dpng -r600 ex4_4.png; system('convert -trim ex4_4.png eps2:ex4_4.eps');
% figure(2), print -dpng -r600 ex4_5.png; system('convert -trim ex4_5.png eps2:ex4_5.eps');
% figure(3), print -dpng -r600 ex4_6.png; system('convert -trim ex4_6.png eps2:ex4_6.eps');

