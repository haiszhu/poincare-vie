function [QuadI,QuadR] = aLapSLVPQuad(p,sI,sR,tbd,tin)
% take in volume mesh irregular boxes (sIr) and regular boxes (sR), 
% evaluate volume potential correction ... 
%   at boundary target (tbd.x), and target within (tin.x).
%
% Hai 06/13/21

tbd.x = tbd.x(:); tin.x = tin.x(:);
rbd = [real(tbd.x(:)),imag(tbd.x(:))]'; rin = [real(tin.x(:)),imag(tin.x(:))]';

%--------------irregular box-------------- 
if numel(sI)>0
    [powY,powX]=meshgrid(0:p-1);
    NbdspaIr = 0; NinspaIr = 0; % for sparse allocation
    for j=1:numel(sI)
        if ~isfield(sI{j},'xToP')
            gridsize = abs(sI{j}.node0(1)-sI{j}.node0(2));
            sc = sI{j}.sc;
            Ztmp = sI{j}.Z;
            nside = numel(sI{j}.node);
            sIrk = sI{j}; sIrk.p = p; sIrk.Z = @(t) 2*(Ztmp(t)-sc)/gridsize;
            sIrk.node = sIrk.Z(linspace(0,2*pi,nside+1)); sIrk.node = sIrk.node(1:end-1);
            xToP = xToPara(sIrk); 
            sI{j}.xToP = xToP; % sIxToP = cellfun(@(d)d.xToP,sI,'uniformoutput',0); save('sIxToP3.mat','sIxToP')
        end
    end
    for j=1:numel(sI)
        gridsize = abs(sI{j}.node0(1)-sI{j}.node0(2));
        sc0 = sI{j}.sc0; % center of this patch
        sc = sI{j}.sc;
        Ztmp = sI{j}.Z;

        % distance criteria (not ideal...)
        idx = abs(tbd.x-sc0) < 1.5*sqrt(2)*gridsize/2;
        tidx = abs(tin.x-sc0) < 1.5*sqrt(2)*gridsize/2;
        xToP = sI{j}.xToP; 

        sI{j}.p = p; %sI{j}.sc0 = sI{j}.xc;

        % close targets...
        tj.x = rbd(1,idx)+1i*rbd(2,idx); sI{j}.tcusp = xToP(2*(tj.x-sc)/gridsize);

        [Utr,sx,sws,coeffFun] = aLapSLVPmatIr(tj,sI{j});
    
        XpowerValf = bsxfun(@power,2*(real(sx)'-real(sc))/gridsize,(0:p-1)'); YpowerValf = bsxfun(@power,2*(imag(sx)'-imag(sc))/gridsize,(0:p-1)'); 
        muTxy = XpowerValf(powX(:)'+1,:).*YpowerValf(powY(:)'+1,:);
        sI{j}.Abdcrr = Utr + 1/(2*pi)*log(abs(bsxfun(@minus,tj.x(:),sx.')))*bsxfun(@times,muTxy',sws(:));
        sI{j}.bdIdx = idx;
        
        tj.x = rin(1,tidx)+1i*rin(2,tidx); sI{j}.tcusp = xToP(2*(tj.x-sc)/gridsize);
        uUtr = aLapSLVPmatIr(tj,sI{j});
        sI{j}.Aincrr = uUtr + 1/(2*pi)*log(abs(bsxfun(@minus,tj.x(:),sx.')))*bsxfun(@times,muTxy',sws(:));
        sI{j}.inIdx = tidx;
        
        sI{j}.sx = sx; sI{j}.sws = sws; sI{j}.coeffFun = coeffFun;
        NbdspaIr = NbdspaIr+numel(sI{j}.Abdcrr);
        NinspaIr = NinspaIr+numel(sI{j}.Aincrr);
    end
    sxIr = cellfun(@(d)d.sx,sI,'uniformoutput',0); % sxIr = vertcat(sxIr{:}); 
    swsIr = cellfun(@(d)d.sws,sI,'uniformoutput',0); swsIr = vertcat(swsIr{:});
    coeffFunIr = cellfun(@(d)d.coeffFun,sI,'uniformoutput',0); 
    % correction matrix from irregular boxes for targets on the boundary
    if numel(tbd.x)>0
        AbdcrrIr = spalloc(numel(tbd.x),numel(sI)*p^2,NbdspaIr); 
        for j=1:numel(sI)
            idx = (j-1)*p^2+(1:p^2);
            AbdcrrIr(sI{j}.bdIdx,idx) = sI{j}.Abdcrr;
        end
    else
        AbdcrrIr = [];
    end
    % correction matrix from irregular boxes for targets within the domain
    if numel(tin.x)>0
        AincrrIr = spalloc(numel(tin.x),numel(sI)*p^2,NinspaIr); 
        for j=1:numel(sI)
            idx = (j-1)*p^2+(1:p^2);
            AincrrIr(sI{j}.inIdx,idx) = sI{j}.Aincrr;
        end
    else
        AincrrIr = [];
    end
else
    AbdcrrIr = []; AincrrIr = []; sxIr = []; coeffFunIr = [];
end

%--------------regular box-------------- 
if numel(sR)>0
    % upsampling...
    [powY,powX]=meshgrid(0:p-1); powX = powX(:)'; powY = powY(:)'; 
    ref = 1; [x0f,w0f]=cheby(ref*p); [XXf,YYf] = meshgrid(x0f);
    % build chebyshev interpolation matrix
    T = 2*cos(pi*((1:p)-1)'*(2*(p:-1:1)-1)/(2*p))/p; T(1,:) = 1/2*T(1,:); % chebyshev approximation
    Tfull = zeros(p^2,p^2);
    for l=1:p
        Tfull((l-1)*p+(1:p),:) = kron(eye(p),T(l,:));
    end
    Tfull = kron(eye(p),T)*Tfull;
    mufTfull = (cos(acos(XXf(:))*powX).*cos(acos(YYf(:))*powY))*Tfull;
    NbdspaR = 0; NinspaR = 0; % for sparse allocation
    for j=1:numel(sR)
        gridsize = abs(sR{j}.node(1)-sR{j}.node(2));
        sc0 = sR{j}.sc; % center of this patch

        % build density approximation on this patch
        sxj = gridsize*(XXf(:)+1i*YYf(:))/2+sc0;
        tws = w0f(:)*w0f(:)'; tws = (gridsize/2)^2*tws(:)';

        % distance criteria
        idx = abs(tbd.x-sc0) < 1.5*sqrt(2)*gridsize/2;
        tidx = abs(tin.x-sc0) < 1.5*sqrt(2)*gridsize/2;
        
        box = sR{j}; box.p = p; box.sc = sR{j}.sc;
        % tbd
        tj.x = tbd.x(idx); %tj.x = tj.x(:);
        [Utr,sx,sws] = aLapSLVPmatR(tj,box);
        sR{j}.Abdcrr = Utr + 1/(2*pi)*bsxfun(@times,log(abs(bsxfun(@minus,tj.x(:),sxj.'))),tws)*mufTfull;
        sR{j}.bdIdx = idx;
        % tin
        tj.x = tin.x(tidx); %tj.x = tj.x(:);
        [uUtr,~,~] = aLapSLVPmatR(tj,box);
        sR{j}.Aincrr = uUtr + 1/(2*pi)*bsxfun(@times,log(abs(bsxfun(@minus,tj.x(:),sxj.'))),tws)*mufTfull;
        sR{j}.inIdx = tidx;
        
        sR{j}.sx = sx; sR{j}.sws = sws;
        NbdspaR = NbdspaR+numel(sR{j}.Abdcrr);
        NinspaR = NinspaR+numel(sR{j}.Aincrr);
    end

    sxR = cellfun(@(d)d.sx,sR,'uniformoutput',0); sxR = vertcat(sxR{:}); 
    swsR = cellfun(@(d)d.sws,sR,'uniformoutput',0); swsR = vertcat(swsR{:}); 
    % correction matrix from regular boxes for targets on the boundary
    if numel(tbd.x)>0
        AbdcrrR = spalloc(numel(tbd.x),numel(sxR),NbdspaR); n0 = 0;
        for j=1:numel(sR)
            idx = n0+(1:numel(sR{j}.sx));
            AbdcrrR(sR{j}.bdIdx,idx) = sR{j}.Abdcrr;
            n0 = n0+numel(sR{j}.sx);
        end
    else
        AbdcrrR = [];
    end
    % correction matrix from regular boxes for targets within the domain
    if numel(tin.x)>0
        AincrrR = spalloc(numel(tin.x),numel(sxR),NinspaR); n0 = 0;
        for j=1:numel(sR)
            idx = n0+(1:numel(sR{j}.sx));
            AincrrR(sR{j}.inIdx,idx) = sR{j}.Aincrr;
            n0 = n0+numel(sR{j}.sx);
        end
    else
        AincrrR = [];
    end
else
    sxR = []; swsR = []; AbdcrrR = []; AincrrR = [];
end

QuadI.Abd = AbdcrrIr; QuadI.Ain = AincrrIr; QuadI.sx = sxIr; QuadI.coeffFun = coeffFunIr; 
QuadR.Abd = AbdcrrR; QuadR.Ain = AincrrR; QuadR.sx = sxR; QuadR.sws = swsR;

end