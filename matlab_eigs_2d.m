function [V,eval] = matlab_eigs_2d(X,K)


[m,n] = size(X);
D = zeros(m,m);
for i = 1:m
    for j = 1:m
        D(i,j) = norm(X(i,:)-X(j,:))^2+eps;
    end
end


scale = 0.04;

WW = exp(-D./scale);

DD=sum(WW,2);

D12=diag(1./sqrt(DD));


% eigs options
    opts.issym = 1;
    opts.isreal = 1;
    opts.tol = 1e-6;
    opts.disp = 0;

 [V,eval] = eigs(D12*WW*D12,K,'lm',opts); 
 
end