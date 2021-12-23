function amat=koorn_coefs2vals(nmax, npols, uvs, polsKN)
%   !
%   ! This routine returns a square matrix that maps coefficients in an
%   ! orthonormal polynomial expansion on the simplex (0,0), (1,0),
%   ! (0,1) to function values at the points uvs. The point locations
%   ! can be arbitrary, but note that they WILL affect the conditioning
%   ! of this matrix.
%   !
%   ! Input:
%   !   nmax, npols - order of expansion, it should be the case that
%   !       npols = (nmax+1)(nmax+2)/2
%   !   uvs - point locations, npols of them
%   !
%   ! Output:
%   !   amat - matrix such that vals = amat*coefs
%   !
%   !
[a,b]=size(uvs);
  amat = zeros(b,(nmax+1)*(nmax+2)/2);
  for i=1:b
      [npols2, pols]=koorn_pols(uvs(:,i), nmax, polsKN);
      amat(i,:) = pols;
  end
end
  


function [npols, pols]=koorn_pols(uv, nmax, polsKN)

%   !
%   ! This subroutine evalutes a bunch of orthogonal polynomials on the
%   ! simplex with vertices (0,0), (1,0), (0,1). The polynomials
%   ! computed by this routine are normalized and rotated Koornwinder
%   ! polynomials, which are classically defined on the triangle with
%   ! vertices (0,0), (1,0), (1,1), and given analytically as:
%   !
%   !   K_{n,k}(x,y) = P_{n-k}^(0,2k+1) (2x+1)  x^k  P_k(2y/x-1)
%   !
%   ! After mapping to the uv-simplex via the change of variables
%   !
%   !      u = y,  v = 1-x
%   !
%   ! these polynomials are given as
%   !
%   !   K_{n,k}(u,v) = P_{n-k}^(0,2k+1) (1-2v)  (1-v)^k  \cdot
%   !        P_k((2u+v-1)/(1-v))
%   !
%   ! See Koornwinder 1975 for details, or the NIST handbook.
%   !
%   ! Input:
%   !   uv - a point in the simplex (0,0), (1,0), (0,1)
%   !   nmax - maximum degree polynomials to compute, a total of
%   !       (nmax+1)(nmax+2)/2 values are returned
%   !
%   ! Output:
%   !   npols - number of pols = (nmax+1)*(nmax+2)/2
%   !   pols - values of all the polynomials, ordered in the following
%   !       (n,k) manner:
%   !
%   !            (0,0), (1,0), (1,1), (2,0), (2,1), (2,2), etc.
%   !
%   !
% 
%  implicit real *8 (a-h,o-z)
%  real *8 :: uv(2), pols(*)

%  real *8 :: legpols(0:100), jacpols(0:100,0:100)
% 
% https://github.com/fastalgorithms/smooth-surface, this is a slightly
% modified version, if something goes wrong, it is for sure my fault. Hai

legpols=zeros(1,nmax+1);
jacpols=zeros(nmax+1,nmax+1);

  done = 1;
  


  u = uv(1);
  v = uv(2);
  z = 2*u+v-1;
  y = 1-v;

  legpols(1+0) = 1;
  legpols(1+1) = z;

  for k=1:nmax-1
    legpols(1+k+1) = ((2*k+1)*z*legpols(1+k) - k*legpols(1+k-1)*y*y)/(k+1);
  end

  x = 1-2*v;
  
  K = 0:nmax;
  beta = 2*K+1;
  jacpols(1+0,1+K) = 1;
  jacpols(1+1,1+K) = (-beta + (2+beta)*x)/2;
  for n = 1:nmax-1
    an = (2*n+beta+1).*(2*n+beta+2)/2/(n+1)./(n+beta+1);
    bn = (-beta.^2).*(2*n+beta+1)/2/(n+1)./(n+beta+1)./(2*n+beta);
    cn = n*(n+beta).*(2*n+beta+2)/(n+1)./(n+beta+1)./(2*n+beta);
    jacpols(1+n+1,:) = (an*x + bn).*jacpols(1+n,:) - cn.*jacpols(1+n-1,:);
  end
  
%   n_idx = 0:nmax; k_idx = 0:nmax;
%   [ntmp,ktmp] = meshgrid(n_idx,k_idx); uptri_idx = logical(triu(ones(nmax+1)));
%   K = ktmp(uptri_idx); N = ntmp(uptri_idx);
  K = polsKN(:,1); N = polsKN(:,2);
  SC = sqrt(done./(2*K+1)./(2*N+2));
  pols = (legpols(1+K)'.*jacpols((1+N-K)+(nmax+1)*K)./SC)';

  npols = numel(pols);
  
end 
