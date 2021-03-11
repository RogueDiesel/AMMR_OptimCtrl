Nc = 4; Nf = 20;
figure
plot(Phi'*y(1:2*Nf))  %x2(flap speed) time plot
hold on
plot(Phi'*y(2*Nf+1:4*Nf) + y(4*Nf+2*Nc+1))  %x3(generator speed) time plot
figure
plot(Phi'*Dphi*y(2*Nf+1:4*Nf)) %dx3(generator acceleration) time plot 
figure
plot(y(4*Nf+1:4*Nf+Nc)) %up plot
hold on
plot(y(4*Nf+Nc+1:4*Nf+2*Nc)); %uq plot, up+uq = abs(u2), up-uq = u2; 
% plot(Te)
