function u=LapSLPCloseApplyp(NPt,s,t,f,side,Sc)
%
%
% Hai 01/21/21

len = s.p;
if isfield(NPt,'s')
if nargin == 5
    u.s = zeros(length(s.x),1); ls = 0;
    for k = 1:s.np
        idx1=([k;NPt.s{k}]-1)*len;idx2=1:len;
        tidx = reshape((idx1*ones(1,len)+ones(numel(idx1),1)*idx2)',[],1);
        ff=f(ls+1:ls+len);
        ss.np=1; ss.p=s.p; ss.x=s.x(ls+1:ls+len); ss.xp=s.xp(ls+1:ls+len); ss.xpp=s.xpp(ls+1:ls+len); ss.xlo=s.xlo(k); ss.xhi=s.xhi(k); ss.t=s.t(ls+1:ls+len); ss.tlo=s.tlo(k); ss.thi=s.thi(k); ss.sp=s.sp(ls+1:ls+len); ss.w=s.w(ls+1:ls+len); ss.nx=s.nx(ls+1:ls+len); ss.ws=s.ws(ls+1:ls+len); ss.wxp=s.wxp(ls+1:ls+len);
        tx.x = s.x(tidx); skip = numel(tx.x);
        Aspecial = LapSLPSpecialMat(tx,ss,side); uSpecial = Aspecial*ff; Anaive = LapSLPmat(tx,ss); Anaive(1:skip+1:end) = 0; uNaive = Anaive*ff; 
        uCrr = uSpecial-uNaive;
        u.s(tidx) = u.s(tidx)+uCrr;
        ls = ls+len;
    end
else
    u.s = Sc.ssparse*f;
end
end

if isfield(NPt,'t')
if nargin == 5
    u.t = zeros(length(t.x),1); ls = 0;
    for k = 1:s.np
        if ~isempty(NPt.t{k})
            ff=f(ls+1:ls+len);
            ss.np=1; ss.p=s.p; ss.x=s.x(ls+1:ls+len); ss.xp=s.xp(ls+1:ls+len); ss.xpp=s.xpp(ls+1:ls+len); ss.xlo=s.xlo(k); ss.xhi=s.xhi(k); ss.t=s.t(ls+1:ls+len); ss.tlo=s.tlo(k); ss.thi=s.thi(k); ss.sp=s.sp(ls+1:ls+len); ss.w=s.w(ls+1:ls+len); ss.nx=s.nx(ls+1:ls+len); ss.ws=s.ws(ls+1:ls+len); ss.wxp=s.wxp(ls+1:ls+len);
            tx.x = t.x(NPt.t{k}); 
            Aspecial = LapSLPSpecialMat(tx,ss,side); uSpecial = Aspecial*ff; Anaive = LapSLPmat(tx,ss); uNaive = Anaive*ff; 
            uCrr = uSpecial-uNaive;
            u.t(NPt.t{k}) = u.t(NPt.t{k})+uCrr;
        end
        ls = ls+len;
    end
else
    u.t = Sc.tsparse*f;
end
end

% keyboard

end