function val=koorn_evalexp(nmax, npols, uv, coefs)

%   !
%   ! Evaluate the orthgonal polynomial expansion with given
%   ! coefficients at the points uv
%   !
%   ! Input:
%   !   nmax, npols - number of terms in expansion,
%   !       npols = (nmax+1)*(nmax+2)/2
%   !   uv - point in the simplex at which to evaluate
%   !   coefs - coefs of expansion, ordered the same as koorn_pols
%   !
%   ! Output:
%   !   val - the value of the expansion at uv
%   !
%   !

  [npols2, pols]=koorn_pols(uv, nmax);
  val = 0;
  for i = 1:npols
    val = val + coefs(i)*pols(i);
  end

end 