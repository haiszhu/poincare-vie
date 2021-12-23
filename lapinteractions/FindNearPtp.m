function NPt=FindNearPtp(s,if_s,t,if_t,dd)

% hacked Barnett version which replaces Stats Toolbox kdtree
% with an open-source version;
% but I don't really understand Gary's code here. 11/30/17
%
% panel version
% s for source, t for target, dd is for example 5 in (5h rule)
% output is based on smallest structure, panel
% NPt.s and NPt.t each has num of panel attributes.
% NPt.s{k} (for kth panel), stores panel indexs close to kth panel (exclude itself)
% NPt.t{k} (for kth panel), stores target point indexs close to kth panel.
% 05/07/18 Hai

Y=[real(s.x),imag(s.x)];
M = s.np;  
len = s.p;

if if_s==1       % self search over s
    NPt.s=cell(1,M);
    tree = kdtree_build(Y);   % Y is Mx2
    ls=0;
    for k=1:M
      idy = [];
      j = (k-1)*len+(1:len); 
      h=max(s.ws(j)); d=dd*h;
      for m=ls+1:ls+len   % can only query one at a time
        idy = [idy; kdtree_ball_query(tree, Y(m,:),d)];  % ahb
      end
      NPt.s{k}=unique(floor((idy-1)/len)+1);
      NPt.s{k}=NPt.s{k}(NPt.s{k}~=k);
      ls=ls+len;
    end
end

if if_t==1      % tree on targs, search each source
    NPt.t=cell(1,M);
    X=[real(t.x),imag(t.x)];
    tree = kdtree_build(X);
    ls=0;
    for k=1:M
        idx = [];
        j = (k-1)*len+(1:len); 
        h=max(s.ws(j)); d=dd*h;
        for m=ls+1:ls+len   % can only query one at a time
          idx = [idx; kdtree_ball_query(tree, Y(m,:),d)];  % ahb
        end
        NPt.t{k}=unique(idx);   % NPt.t{k} corresponding to targets needs to be corrected for kth panel
        ls=ls+len;
    end
end
