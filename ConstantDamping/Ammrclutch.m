% The below parameters works best sofar (This case does converge and make some sense)
T = 10;     Nc = 40;      Nf = 20;      Ff = 1e3; 

tk = linspace(0, T, Nc+2);
tk = tk(2:end-1);
Te = 1000*sin(pi/5*tk); 


w0 = 2*pi/T;
W = w0*(1:Nf)';
Wtk = W*tk;
Phi = zeros(2*size(Wtk,1),size(Wtk,2));             
Phi(1:2:end,:) = cos(Wtk);
Phi(2:2:end,:) = sin(Wtk);

% Building derivative matrix
d = zeros(2*Nf-1,1);
d(1:2:end) = 1:Nf;
Dphi = diag(d,1);
Dphi = w0*(Dphi - Dphi');
Iphi = inv(Dphi);

R = @(x) fourierclutch(x,Phi',Dphi,Iphi,Te',Nc,Nf,Ff);
Fobj = @(x) netpower(x,Nc,Nf);

opts = optimoptions(@fmincon,'Algorithm','sqp','Diagnostics','on','Display','iter',...
'MaxFunctionEvaluations',inf,'MaxIterations',inf,'OutputFcn',[],...
'PlotFcn',{'optimplotfval','optimplotconstrviolation','optimplotstepsize','optimplotfirstorderopt'},...
'TolFun',1e-12,'TolCon',1e-12,'TolX',1e-14,'GradObj','on','GradConstr','on','CheckGradients',true,...
'FiniteDifferenceType','central');

%Initialization
x0 = zeros(4*Nf+2*Nc+1,1);
x0(1:4*Nf) = -10 + 20*rand(4*Nf,1);  %Fourier coefficients
x0(4*Nf+1:4*Nf+2*Nc) = rand(2*Nc,1); %Slack variable up,uq
x0(4*Nf+2*Nc+1)= 10*rand;  % x3 dc term

%Bounds
lb = zeros(4*Nf+2*Nc+1,1); ub = zeros(4*Nf+2*Nc+1,1);

lb(1:4*Nf)=-10; ub(1:4*Nf)=10;     % x2 x3
lb(4*Nf+1:4*Nf+2*Nc)= 0; ub(4*Nf+1:4*Nf+2*Nc)= 1;   %slack variables
lb(4*Nf+2*Nc+1)= 0; ub(4*Nf+2*Nc+1)= 10; %x3 dc term

%Linear constraint keeping x3 nonnegative
A = zeros(Nc,4*Nf+2*Nc+1);
A(:,2*Nf+1:4*Nf) = Phi';      A(:,4*Nf+2*Nc+1) = 1;
A = -A; b = zeros(Nc,1);
[y, fval] = fmincon(Fobj,x0,A,b,[],[],lb,ub,R,opts);
