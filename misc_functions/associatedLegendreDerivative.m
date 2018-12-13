function [dP] = associatedLegendreDerivative(n,x)
%ASSOCIATEDLEGENDREDERIVATIVE
% First derivative of associated legendre function using recursion formula
% and boundary conditions
%

m = -n:n;
M = repmat(m,length(x),1);
X = repmat(x,1,length(m));
Pnm = associatedLegendre(n,x);

DP1 = Pnm;

% Handle m=0
i = find(m==0);
DP2 = [Pnm(:,i)./sqrt(1-x.^2) Pnm(:,i+1:end)];
DP2 = [fliplr(DP2(:,2:end)) DP2];

dP = 1./(1-X.^2) .* (n*x.*DP1 -(n+M).*(n-M+1).*sqrt(1-X.^2).*DP2);

% Handle boundaries at x = +/-1
i1 = find(abs(x+1)<0.0000001);
i2 = find(abs(x-1)<0.0000001);

k0 = find(abs(m)==0);
k1 = find(abs(m)==1);
k2 = find(abs(m)==2);
k3 = find(abs(m)>=3);

dP(i1,k0) = n*(n+1)/2;
dP(i1,k1) = Inf;
dP(i1,k2) = -(n-1)*n*(n+1)*(n+2)/4;
dP(i1,k3) = 0;

dP(i2,k0) = (-1)^(n+1)*n*(n+1)/2;
dP(i2,k1) = (-1)^n*Inf;
dP(i2,k2) = (-1)^n*(n-1)*n*(n+1)*(n+2)/4;
dP(i2,k3) = 0;

end

