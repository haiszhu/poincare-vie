function u=As_operator(X,s,NPt,side,Sc)
% define operator...
% X is vector in GMRES
%
% Hai 01/21/21

% close interaction...
NPts={}; NPts.s=NPt.s;
uSclose=LapSLPCloseApplyp(NPts,s,s,X,side,Sc);
% far interaction
uSfar=LapSLPFMM2D(s,1,s,0,X,4);
u=uSclose.s(:)+uSfar(:);
% keyboard

end