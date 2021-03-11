function [f, grad] = netpower(x,Nc,Nf)
v = x(2*Nf+1:4*Nf); vmean = x(4*Nf+2*Nc+1); 
% u2p = x(4,1:end-1); u2q = x(5,1:end-1);
f = -5*50*sum(v.^2) - 10*50*vmean^2;
% f = sum(u2p.*u2q);
if nargout > 1
    grad  = x;
    grad(1:2*Nf) = 0; grad(2*Nf+1:4*Nf) = -500*v; grad(4*Nf+1:4*Nf+2*Nc) = 0;
    grad(4*Nf+2*Nc+1) = -1000*vmean;
end