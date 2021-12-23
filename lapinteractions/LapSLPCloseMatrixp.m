function [Sc,Ss]=LapSLPCloseMatrixp(NPt,s,t,side)
% Sc is for self or target correction matrix... Ss is self for block
% diagonal preconditioner...
%
%
% Hai 01/21/21

len = s.p; M = s.np;
if isfield(NPt,'s')
    Sc.s=cell(1,M);Sc.i=cell(1,M);Sc.j=cell(1,M); ls = 0;
    if nargout>1, Sself.s=cell(1,M);Sself.i=cell(1,M);Sself.j=cell(1,M); end
    for k = 1:M
        idx1=([k;NPt.s{k}]-1)*len;idx2=1:len;
        tidx = reshape((idx1*ones(1,len)+ones(numel(idx1),1)*idx2)',[],1);
        ss.np=1; ss.p=s.p; ss.x=s.x(ls+1:ls+len); ss.xp=s.xp(ls+1:ls+len); ss.xpp=s.xpp(ls+1:ls+len); ss.xlo=s.xlo(k); ss.xhi=s.xhi(k); ss.t=s.t(ls+1:ls+len); ss.tlo=s.tlo(k); ss.thi=s.thi(k); ss.sp=s.sp(ls+1:ls+len); ss.w=s.w(ls+1:ls+len); ss.nx=s.nx(ls+1:ls+len); ss.ws=s.ws(ls+1:ls+len); ss.wxp=s.wxp(ls+1:ls+len);
        tx.x = s.x(tidx); skip = numel(tx.x);
        Aspecial = LapSLPSpecialMat(tx,ss,side); Anaive = LapSLPmat(tx,ss); Anaive(1:skip+1:end) = 0; 
        Sc.s{k} = Aspecial - Anaive;
        rows = tidx*ones(1,len); columns = ones(size(tidx))*[ls+1:ls+len];
        Sc.i{k} = rows(:); Sc.j{k} = columns(:); Sc.s{k} = Sc.s{k}(:);
        if nargout>1
            Aspecial_s = LapSLPSpecialMat(ss,ss,side);
            rows_s=[ls+1:ls+len]'*ones(1,len); columns_s=ones(len,1)*[ls+1:ls+len]; 
            Sself.i{k}=rows_s(:); Sself.j{k}=columns_s(:); Sself.s{k}=Aspecial_s(:);
        end
        ls = ls+len;
    end
    Sci=cell2mat(Sc.i');Scj=cell2mat(Sc.j');Scs=cell2mat(Sc.s');
    Sc = rmfield(Sc,{'i','j','s'});
    Sc.ssparse=sparse(Sci,Scj,Scs,M*len,M*len);
    if nargout>1
        Ssi=cell2mat(Sself.i');Ssj=cell2mat(Sself.j');Sss=cell2mat(Sself.s');
        Ss=sparse(Ssi,Ssj,Sss,M*len,M*len);
    end
end

if isfield(NPt,'t')
    Sc.t=cell(1,M);Sc.i=cell(1,M);Sc.j=cell(1,M); ls = 0;
    for k = 1:M
        if ~isempty(NPt.t{k})
            ss.np=1; ss.p=s.p; ss.x=s.x(ls+1:ls+len); ss.xp=s.xp(ls+1:ls+len); ss.xpp=s.xpp(ls+1:ls+len); ss.xlo=s.xlo(k); ss.xhi=s.xhi(k); ss.t=s.t(ls+1:ls+len); ss.tlo=s.tlo(k); ss.thi=s.thi(k); ss.sp=s.sp(ls+1:ls+len); ss.w=s.w(ls+1:ls+len); ss.nx=s.nx(ls+1:ls+len); ss.ws=s.ws(ls+1:ls+len); ss.wxp=s.wxp(ls+1:ls+len);
            tx.x = t.x(NPt.t{k}); 
            Aspecial = LapSLPSpecialMat(tx,ss,side); Anaive = LapSLPmat(tx,ss); 
            Sc.t{k} = Aspecial - Anaive;
            rows = NPt.t{k}*ones(1,len); columns = ones(size(NPt.t{k}))*[ls+1:ls+len];
            Sc.i{k} = rows(:); Sc.j{k} = columns(:); Sc.t{k}=Sc.t{k}(:);
        end
        ls = ls+len;
    end
    Sci=cell2mat(Sc.i');Scj=cell2mat(Sc.j');Sct=cell2mat(Sc.t');
    Sc = rmfield(Sc,{'i','j','t'});
    Sc.tsparse=sparse(Sci,Scj,Sct,numel(t.x),M*len);
end


% keyboard

end