% <============================================================================>
% < Author: Gerardo Chowell  ==================================================>
% <============================================================================>

function dx=GRM(t,x,params0,extra0)

% Generalized Richards model (GrowthPredict flag1=2)
% parameters in order: r, p, a, K

dx=zeros(1,1);

dx(1,1)=params0(1)*(x(1,1).^params0(2)).*(1-(x(1,1)/params0(4)).^params0(3));
