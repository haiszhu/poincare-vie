function [A,T] = StoSLPmat(t,s) % single-layer kernel matrix
% t = target seg (x,nx cols), s = src seg
% BW 4/9/19
N = numel(s.x); M = numel(t.x);
d = repmat(t.x, [1 N]) - repmat(s.x.', [M 1]);    % C-# displacements mat
r = abs(d);    % dist matrix R^{MxN}

log_part = -log(r);

d1 = real(d)./r;
d2 = imag(d)./r;
if numel(s.x)==numel(t.x) && max(abs(s.x-t.x))<1e-14
    d1(diagind(d1)) = real(s.tang);
    d2(diagind(d2)) = imag(s.tang);
end  % self? diagonal term for the cross part

cross12 = d1.*d2; % off-diag block of the cross part
A = [log_part + d1.^2, cross12; cross12, log_part + d2.^2];
A = bsxfun(@times, A, repmat(s.ws(:).'/4/pi, [1 2]));
if nargout > 1
    nx = repmat(t.nx, [1 N]);
    dot_part = -real(nx./d);
    if numel(s.x)==numel(t.x) && max(abs(s.x-t.x))<1e-14
        dot_part(diagind(dot_part)) = -s.cur/2;
    end     % self? correct diagonal limits
    cross12 = dot_part.*d1.*d2; % off-diag block
    T = [dot_part.*d1.^2, cross12; cross12, dot_part.*d2.^2];
    T = bsxfun(@times, T, repmat(s.ws(:).'/pi, [1 2]));
end
end
