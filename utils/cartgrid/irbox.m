function [sIr, sIrIdxMat] = irbox(xgrid,ygrid,s,tIC,side,tol)
% assuming gridsize is the same, square boxes...
% one at a time
%
% Hai 05/03/21, modified 05/23/21 to return index matrix of size (numel(ygrid)-1) X (numel(xgrid)-1)

sIrIdxMat = NaN(numel(ygrid)-1,numel(xgrid)-1); xgridc0 = 1/2*(xgrid(1)+xgrid(2)); ygridc0 = 1/2*(ygrid(1)+ygrid(2)); % box centers anchored here

gridsize = xgrid(2)-xgrid(1);
if nargin < 5, side = 'i'; end
if nargin < 6, tol = 1e-10; end
pflag = ischar(s);  % take string input from python...
if ~isfield(s,'Z') 
    if pflag
        s = str2func(s);
    end
    s = struct('Z',s); 
end
delta = 0.00001; tmp = s.Z(linspace(0,2*pi,1000)); tmpsc = mean(tmp);
[Xgrid,Ygrid] = meshgrid(xgrid,ygrid); 
if side == 'i'
    [IN,~] = inpolygon(Xgrid,Ygrid,(1+delta)*real(tmp-tmpsc)+real(tmpsc),(1+delta)*imag(tmp-tmpsc)+imag(tmpsc)); % using delta might not be a good idea
    INgrid = IN(1:end-1,1:end-1)&IN(2:end,1:end-1)&IN(1:end-1,2:end)&IN(2:end,2:end);
    ONgrid = (IN(1:end-1,1:end-1)|IN(2:end,1:end-1)|IN(1:end-1,2:end)|IN(2:end,2:end))&~INgrid; % at least one corner within s, but not all four corners
else
    [IN,~] = inpolygon(Xgrid,Ygrid,(1-delta)*real(tmp-tmpsc)+real(tmpsc),(1-delta)*imag(tmp-tmpsc)+imag(tmpsc));
    nINgrid = ~IN(1:end-1,1:end-1)&~IN(2:end,1:end-1)&~IN(1:end-1,2:end)&~IN(2:end,2:end);
    ONgrid = (~IN(1:end-1,1:end-1)|~IN(2:end,1:end-1)|~IN(1:end-1,2:end)|~IN(2:end,2:end)) & ~nINgrid; % at least one corner outside s
end
[onIdy,onIdx] = find(ONgrid==1);

sIr = {}; sxIC = s.Z(tIC); j = 0;
for k=1:numel(onIdx)
    sIrksc = 1/2*(xgrid(onIdx(k))+xgrid(onIdx(k)+1)) + 1i/2*(ygrid(onIdy(k))+ygrid(onIdy(k)+1));
    sIrkZ0 =@(t) sIrksc + (((t>=0)&(t<pi/2)).*(1./(sin(t)+cos(t))).*(cos(t)+1i*sin(t))*exp(1i*pi/4)...
                     + ((t>=pi/2)&(t<pi)).*(1./(sin(t-pi/2)+cos(t-pi/2))).*(cos(t-pi/2)+1i*sin(t-pi/2))*exp(1i*3*pi/4)...
                     + ((t>=pi)&(t<3*pi/2)).*(1./(sin(t-pi)+cos(t-pi))).*(cos(t-pi)+1i*sin(t-pi))*exp(1i*5*pi/4)...
                     + ((t>=3*pi/2)&(t<=2*pi)).*(1./(sin(t-3*pi/2)+cos(t-3*pi/2))).*(cos(t-3*pi/2)+1i*sin(t-3*pi/2))*exp(1i*7*pi/4))*gridsize/sqrt(2);
    sIrknode0 = sIrkZ0(linspace(0,2*pi,5)); % this is not really needed after setting up irregular box, neither is sIr{k}.Z0 
    sIrknode0 = sIrknode0(1:4);
    
    % regular box nodes, counterclockwise order
    INk = [IN(onIdy(k)+1,onIdx(k)+1),IN(onIdy(k)+1,onIdx(k)),IN(onIdy(k),onIdx(k)),IN(onIdy(k),onIdx(k)+1)];
    if side == 'e', INk = ~INk; end
    DnodeR = sIrknode0(INk); 
    if prod(INk==[1,0,1,1]), DnodeR = [DnodeR(2:3),DnodeR(1)]; end % correction of order for a few patterns
    if prod(INk==[1,1,0,1]), DnodeR = [DnodeR(3),DnodeR(1:2)]; end
    if prod(INk==[1,0,0,1]), DnodeR = [DnodeR(2),DnodeR(1)]; end
    
    % irregular box nodes, counterclockwise order
    tmp1 = real((sxIC-sIrknode0(1))'*(sIrknode0(2)-sIrknode0(1))); % anchored at northeast sIr{k}.node0(1), check inner product
    tmp2 = real((sxIC-sIrknode0(1))'*(sIrknode0(4)-sIrknode0(1))); % see which interception tIC belongs to this box
    inIC = (tmp1<=norm(sIrknode0(2)-sIrknode0(1))^2+tol) & (tmp1>=-tol) & (tmp2<=norm(sIrknode0(4)-sIrknode0(1))^2+tol) & (tmp2>=-tol);
    if inIC(end)&&inIC(end-1), inIC(1) = 0; end % if tIC = 0, and 2*pi account for two, check 2nd to last one. (tIC was processed by 'uniquetol()', 
    DnodeIr = sxIC(inIC); tkIC = tIC(inIC);     % assuming mesh is reasonable, then sum(inIC)=2 except for the case t=0, 2*pi)
    if numel(DnodeIr)>2, 1; end     % do something else here, if enter and leave the box more than once...
    
    % try to take care of shared node between DnodeR and DnodeIr
    shaNodeFlag = sum(abs(DnodeIr - DnodeR) < tol);
    if sum(shaNodeFlag) > 0
        DnodeR = DnodeR(~shaNodeFlag);
    end
    
    if numel(DnodeIr)==2 % both enter and leave the box 
        if side == 'e', DnodeIr = DnodeIr(2:-1:1); tkIC = tkIC(2:-1:1); end
        j = j+1;
        DnodeR = DnodeR(~sum(abs(DnodeIr-DnodeR)<tol)); % there might be repeated nodes, shared between DnodeIr and DnodeR
        if tkIC(2)-tkIC(1)>pi && side == 'i'
            tkIC = [tkIC(2); 2*pi+tkIC(1)]; DnodeIr = [DnodeIr(2);DnodeIr(1)]; 
        end % some rule for potential box that has part of s curve containing s.Z(0), s.Z(2*pi)
        if tkIC(1)-tkIC(2)>pi && side == 'e'
            tkIC = [2*pi+tkIC(2); tkIC(1)]; DnodeIr = [DnodeIr(2);DnodeIr(1)]; 
        end
        % define boundary
        if numel(DnodeR)==0
            sIr{j}.Z = @(t)   ((t>=0)&(t<pi)).*(DnodeIr(2)+t/pi*(DnodeIr(1)-DnodeIr(2)))...
                            + ((t>=pi)&(t<=2*pi)).*s.Z(tkIC(1)+(t-pi)/pi*(tkIC(2)-tkIC(1)));
            sIr{j}.node = sIr{j}.Z([0,pi]);
        elseif numel(DnodeR)==1
            sIr{j}.Z = @(t)   ((t>=0)&(t<2*pi/3)).*(DnodeIr(2)+t/(2*pi/3)*(DnodeR-DnodeIr(2)))...
                            + ((t>=2*pi/3)&(t<4*pi/3)).*(DnodeR+(t-2*pi/3)/(2*pi/3)*(DnodeIr(1)-DnodeR))...
                            + ((t>=4*pi/3)&(t<=6*pi/3)).*s.Z(tkIC(1)+(t-4*pi/3)/(2*pi/3)*(tkIC(2)-tkIC(1)));
            sIr{j}.node = sIr{j}.Z([0,2*pi/3,4*pi/3]);
        elseif numel(DnodeR)==2
            sIr{j}.Z = @(t)   ((t>=0)&(t<pi/2)).*(DnodeIr(2)+t/(pi/2)*(DnodeR(1)-DnodeIr(2)))...
                            + ((t>=pi/2)&(t<pi)).*(DnodeR(1)+(t-pi/2)/(pi/2)*(DnodeR(2)-DnodeR(1)))...
                            + ((t>=pi)&(t<3*pi/2)).*(DnodeR(2)+(t-pi)/(pi/2)*(DnodeIr(1)-DnodeR(2)))...
                            + ((t>=3*pi/2)&(t<=2*pi)).*s.Z(tkIC(1)+(t-3*pi/2)/(pi/2)*(tkIC(2)-tkIC(1)));
            sIr{j}.node = sIr{j}.Z([0,pi/2,pi,3*pi/2]);
        elseif numel(DnodeR)==3
            sIr{j}.Z = @(t)   ((t>=0)&(t<2*pi/5)).*(DnodeIr(2)+t/(2*pi/5)*(DnodeR(1)-DnodeIr(2)))...
                            + ((t>=2*pi/5)&(t<4*pi/5)).*(DnodeR(1)+(t-2*pi/5)/(2*pi/5)*(DnodeR(2)-DnodeR(1)))...
                            + ((t>=4*pi/5)&(t<6*pi/5)).*(DnodeR(2)+(t-4*pi/5)/(2*pi/5)*(DnodeR(3)-DnodeR(2)))...
                            + ((t>=6*pi/5)&(t<8*pi/5)).*(DnodeR(3)+(t-6*pi/5)/(2*pi/5)*(DnodeIr(1)-DnodeR(3)))...
                            + ((t>=8*pi/5)&(t<=10*pi/5)).*s.Z(tkIC(1)+(t-8*pi/5)/(2*pi/5)*(tkIC(2)-tkIC(1)));
            sIr{j}.node = sIr{j}.Z([0,2*pi/5,4*pi/5,6*pi/5,8*pi/5]);
        else
            disp('      maybe some edge case...')
            disp(' ')
            keyboard
        end
        sIr{j}.tIC = tkIC(1:2)'; sIr{j}.sc = sIrksc; sIr{j}.Z0 = sIrkZ0; sIr{j}.node0 = sIrknode0; 
        scxIdx = round((real(sIr{j}.sc)-xgridc0)/gridsize)+1; scyIdx = round((imag(sIr{j}.sc)-ygridc0)/gridsize)+1; 
        sIrIdxMat(scyIdx,scxIdx) = j;
    end
    
end
if pflag
    node0 = cellfun(@(d)d.node0,sIr,'uniformoutput',0); node0 = vertcat(node0{:}); 
    tkIC = cellfun(@(d)d.tIC,sIr,'uniformoutput',0); tkIC = vertcat(tkIC{:}); 
    node = cellfun(@(d)d.node,sIr,'uniformoutput',0); 
    for k=1:numel(node)
        node{k} = [node{k},nan(1,5-numel(node{k}))];
    end
    node = vertcat(node{:});
    sIr = [node0,tkIC,node]; % of size num of irregular boxes X 9 (1st 4 columns are regular nodes, next 2 columns are the intersection)
end

end
