function [c,ceq,gradc,gradceq] = fourierclutch(x,Phi,Dphi,Iphi,T,Nc,Nf,Ff)
x2 = x(1:2*Nf); x3 = x(2*Nf+1:4*Nf); 
u2p = x(4*Nf+1:4*Nf+Nc); x3mean = x(4*Nf+2*Nc+1);
u2q = x(4*Nf+Nc+1:4*Nf+2*Nc);
% xs = diag(Iphi,-1);
% ing = sum(x2(2:2:end).*xs(1:2:end));

% ceq = [(Phi*Dphi*x2*500 + 1500*Phi*x2 + 2000*(Phi*Iphi*x2) - T - 1e3*(u2.*(Phi*x3)-abs(u2).*(Phi*x2)));
%     (Phi*Dphi*x3*100 - Phi*u1 + 10*Phi*x3 + 1e3*(abs(u2).*(Phi*x3)-u2.*(Phi*x2)))];
c = [-Phi*x2.*u2p; Phi*x2.*u2q]; 
ceq = [(Phi*Dphi*x2*500 + 1500*Phi*x2 + 2000*(Phi*Iphi*x2) - T - Ff*((u2p-u2q).*(Phi*x3+x3mean)-(u2p+u2q).*(Phi*x2)));
    (Phi*Dphi*x3*100 + (300*(Phi*x3 + x3mean)) + Ff*((u2p+u2q).*(Phi*x3+x3mean)-(u2p-u2q).*(Phi*x2)));
    u2p.*u2q];
if nargout > 2
    gradc = [-Phi.*u2p, zeros(Nc,2*Nf), diag(-Phi*x2), zeros(Nc), zeros(Nc,1);
        Phi.*u2q, zeros(Nc,2*Nf), zeros(Nc), diag(Phi*x2), zeros(Nc,1)]';
    
    gradceq = [Phi*Dphi*500 + 1500*Phi + 2000*Phi*Iphi + Ff*(u2p+u2q).*Phi,...
        -Ff*(u2p-u2q).*Phi, diag(-Ff*(Phi*x3+x3mean)+Ff*Phi*x2),...
        diag(Ff*(Phi*x3+x3mean)+Ff*(Phi*x2)), -Ff*(u2p-u2q);
        -Ff*(u2p-u2q).*Phi, Phi*Dphi*100 - 300*Phi + Ff*(u2p+u2q).*Phi, diag(Ff*(Phi*x3+x3mean)-Ff*(Phi*x2)),...
        diag(Ff*(Phi*x3+x3mean)+ Ff*(Phi*x2)), Ff*(u2p+u2q)-300;
        zeros(Nc,4*Nf), diag(u2q), diag(u2p), zeros(Nc,1)]';
end