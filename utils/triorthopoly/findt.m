function t = findt(tx,sc)
% the order is different... not sure which one to use yet
%
% Hai 12/21/21

vert = [[0,1,0];[0,0,1]]'; % maybe modify this part later to be consistent with (xi,eta) space
x1 = vert(1,1); x2 = vert(2,1); x3 = vert(3,1);
y1 = vert(1,2); y2 = vert(2,2); y3 = vert(3,2);
if nargin < 2, sc = 1/3+1i*1/3; end

% 1st 
v1 = x2 + 1i*y2 - sc; v2 = x3 + 1i*y3 - sc;
a1 =  ((real(tx)*imag(v2)-imag(tx)*real(v2)) - (real(sc)*imag(v2)-imag(sc)*real(v2)))./(real(v1)*imag(v2)-imag(v1)*real(v2));
b1 = -((real(tx)*imag(v1)-imag(tx)*real(v1)) - (real(sc)*imag(v1)-imag(sc)*real(v1)))./(real(v1)*imag(v2)-imag(v1)*real(v2));
idx1 = a1>=0 & b1>=0;
tx1 = tx(idx1);
xtmp1 = ((1-imag(tx1)).*(real(tx1)-real(sc)) + real(tx1).*(imag(tx1)-imag(sc)))./(real(tx1)+imag(tx1)-real(sc)-imag(sc));
ytmp1 = (imag(tx1).*(real(tx1)-real(sc)) + (1-real(tx1)).*(imag(tx1)-imag(sc)))./(real(tx1)+imag(tx1)-real(sc)-imag(sc));
ttmp1 = (1-xtmp1)*2*pi/3 + 4*pi/3;

% 2nd
v1 = x3 + 1i*y3 - sc; v2 = x1 + 1i*y1 - sc;
a2 =  ((real(tx)*imag(v2)-imag(tx)*real(v2)) - (real(sc)*imag(v2)-imag(sc)*real(v2)))./(real(v1)*imag(v2)-imag(v1)*real(v2));
b2 = -((real(tx)*imag(v1)-imag(tx)*real(v1)) - (real(sc)*imag(v1)-imag(sc)*real(v1)))./(real(v1)*imag(v2)-imag(v1)*real(v2));
idx2 = a2>=0 & b2>=0;
tx2 = tx(idx2);
ytmp2 = -real(tx2).*(imag(tx2)-imag(sc))./(real(tx2)-real(sc)) + imag(tx2);
ttmp2 = (1-ytmp2)*2*pi/3;

% 3rd
v1 = x1 + 1i*y1 - sc; v2 = x2 + 1i*y2 - sc;
a3 =  ((real(tx)*imag(v2)-imag(tx)*real(v2)) - (real(sc)*imag(v2)-imag(sc)*real(v2)))./(real(v1)*imag(v2)-imag(v1)*real(v2));
b3 = -((real(tx)*imag(v1)-imag(tx)*real(v1)) - (real(sc)*imag(v1)-imag(sc)*real(v1)))./(real(v1)*imag(v2)-imag(v1)*real(v2));
idx3 = a3>=0 & b3>=0;
tx3 = tx(idx3);
xtmp3 = -imag(tx3).*(real(tx3)-real(sc))./(imag(tx3)-imag(sc)) + real(tx3);
ttmp3 = xtmp3*2*pi/3 + 2*pi/3;

t = NaN(size(tx));
t(idx1) = ttmp1; t(idx2) = ttmp2; t(idx3) = ttmp3;

end

