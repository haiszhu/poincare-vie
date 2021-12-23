function [sR, sRIdxMat] = rbox(xgrid,ygrid,s,si)
% assuming gridsize is the same, square boxes...
% assuming interior of s, exterior of si{k}
% do all at once
%
% Hai 05/03/21, modified 05/23/21 to return index matrix of size (numel(ygrid)-1) X (numel(xgrid)-1)

sRIdxMat = NaN(numel(ygrid)-1,numel(xgrid)-1); xgridc0 = 1/2*(xgrid(1)+xgrid(2)); ygridc0 = 1/2*(ygrid(1)+ygrid(2)); % box centers anchored here

gridsize = xgrid(2)-xgrid(1);
if nargin < 4 
    si = []; 
elseif ischar(si)
    eval(append(si,';'))
end
pflag = ischar(s);  % take string input from python...
if ~isfield(s,'Z') 
    if pflag
        s = str2func(s);
    end
    s = struct('Z',s); 
end

delta = 0.0001; tmp = s.Z(linspace(0,2*pi,1000)); tmpsc = mean(tmp);
[Xgrid,Ygrid] = meshgrid(xgrid,ygrid);  
[IN,~] = inpolygon(Xgrid,Ygrid,(1+delta)*real(tmp-tmpsc)+real(tmpsc),(1+delta)*imag(tmp-tmpsc)+imag(tmpsc)); % using delta might not be a good idea
INgrid = IN(1:end-1,1:end-1)&IN(2:end,1:end-1)&IN(1:end-1,2:end)&IN(2:end,2:end);
if numel(si)==0
    [inIdy,inIdx] = find(INgrid==1);
else
    nINigrid = true(numel(ygrid)-1,numel(xgrid)-1);
    for i=1:numel(si)
        tmpi = si{i}.Z(linspace(0,2*pi,1000)); tmpisc = mean(tmpi);
        [INi,~] = inpolygon(Xgrid,Ygrid,(1-delta)*real(tmpi-tmpisc)+real(tmpisc),(1-delta)*imag(tmpi-tmpisc)+imag(tmpisc));
        nINigrid = nINigrid & (~INi(1:end-1,1:end-1)&~INi(2:end,1:end-1)&~INi(1:end-1,2:end)&~INi(2:end,2:end));
    end
    [inIdy,inIdx] = find((nINigrid&INgrid)==1);
end
sR = {};
for k=1:numel(inIdx)
    sc = 1/2*(xgrid(inIdx(k))+xgrid(inIdx(k)+1)) + 1i/2*(ygrid(inIdy(k))+ygrid(inIdy(k)+1));
    sR{k}.sc = sc;
    sR{k}.Z =@(t) sc + (((t>=0)&(t<pi/2)).*(1./(sin(t)+cos(t))).*(cos(t)+1i*sin(t))*exp(1i*pi/4)...
                     + ((t>=pi/2)&(t<pi)).*(1./(sin(t-pi/2)+cos(t-pi/2))).*(cos(t-pi/2)+1i*sin(t-pi/2))*exp(1i*3*pi/4)...
                     + ((t>=pi)&(t<3*pi/2)).*(1./(sin(t-pi)+cos(t-pi))).*(cos(t-pi)+1i*sin(t-pi))*exp(1i*5*pi/4)...
                     + ((t>=3*pi/2)&(t<=2*pi)).*(1./(sin(t-3*pi/2)+cos(t-3*pi/2))).*(cos(t-3*pi/2)+1i*sin(t-3*pi/2))*exp(1i*7*pi/4))*gridsize/sqrt(2);
    sR{k}.node = sR{k}.Z(linspace(0,2*pi,5));
    sR{k}.node = sR{k}.node(1:4);
    scxIdx = round((real(sc)-xgridc0)/gridsize)+1; scyIdx = round((imag(sc)-ygridc0)/gridsize)+1; 
    sRIdxMat(scyIdx,scxIdx) = k;
end
if pflag
    sR = cellfun(@(d)d.node,sR,'uniformoutput',0); sR = vertcat(sR{:}); 
end

end