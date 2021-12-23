function u=Ast_operator(X,s,t,NPt,side,Sc)
% define operator...
% X is vector in GMRES
%
% Hai 01/21/21

% close interaction...
NPtt={}; NPtt.t=NPt.t;
uSclose=LapSLPCloseApplyp(NPtt,s,t,X,side,Sc);
% far interaction
[~,uSfar]=LapSLPFMM2D(s,0,t,1,X,4);
u=uSclose.t(:)+uSfar(:);
% keyboard

end