function [umatr,vmatr]=koorn_vals2coefs_coefs2vals(korder,kpols,xys,polsKN)
%   !
%   !!   compute vals2coefs and coefs2vals matrix 
%   !
%   !    input
%   !    korder     in: integer
%   !                 order of rokhlin vioreanu (rv) nodes at which function is
%   !                 sampled
%   !    kpols      in: integer
%   !                 number of nodes corresponding to korder rv nodes
%   !
%   !    output
%   !    umatr      out: real *8 (kpols,kpols)
%   !               vals2coefs matrix
%   ! 
%   !    vmatr      out: real *8 (kpols,kpols)
%   !               coefs2vals matrix
  if nargin < 3
    [xys,wts]=get_vioreanu_nodes(korder);
    n_idx = 0:korder; k_idx = 0:korder;
    [ntmp,ktmp] = meshgrid(n_idx,k_idx); uptri_idx = logical(triu(ones(korder+1)));
    K = ktmp(uptri_idx); N = ntmp(uptri_idx);
    polsKN = [K,N];
  end
  vmatr=koorn_coefs2vals(korder,kpols,xys,polsKN);
  umatr=inv(vmatr);
end
