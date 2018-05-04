function [V,D,Rank] = main_aca(X,evalnz,MaxRank)
    [w,Rank] = aca_sym2(X,MaxRank);
    N = size(w,2);
    unitvec = ones(N,1);

    D12 = w'*(w*unitvec);
    D12 = 1./sqrt(D12);
    D12 = sparse(1:N,1:N,D12);
    [Q,R] = qr(D12*w',0);
    [V,D] = eigs(R*R',evalnz);
    D = diag(D);
    V = Q*V;
end

function [a,nu] = aca_sym2(u,K)
%% INPUT    Image as n times N matrix in grayscale
%%
%% OUTPUT   nu rank approximation of the weight matrix corresponding to the 
%%          Graph Laplace Operator

% Perform adaptive cross approximations
[n,~]=size(u);
N = n;
M = n;
R = ones(M,1);

% K is number of max rank
%% pick inital point for aca

% ik = randi(M,1,1);
ik = 1;         % use arbitrary initial index


% set tolerance for stop criterion
tol = 1e-6;

for nu=1:K
    %% acaf for graph laplace spectral image clustering
    a(nu,:) = returnsim(u,ik,M);

    %%  substract the rank one approx
    for mu=1:nu-1
        a(nu,:)=a(nu,:)-a(mu,ik)*a(mu,:);
    end
    
    %compute the maximal entry in modulus
    [~,ik]=max(abs(a(nu,:)));
    delta = a(nu,ik);
    %D(nu) = delta;
    
    %fprintf('delta = %f\n',abs(delta));
    
    % the algorithm terminates with the exact rank [nu-1]
    if max(R)<tol
        a = a(1:nu-1,:);
        break;
    end
    if delta<tol
        % fprintf('approx matrix of rank %d\n',nu-1);
        a = a(1:nu-1,:);
        break
    else
        a(nu,:) = 1/sqrt(delta)*a(nu,:);
        for i = 1:M     % we only need the diagonal of R
            %R(i,i) = R(i,i) - a(:,i)'*a(:,i);
            R(i) = R(i) - a(nu,i).^2;
        end
        %R(:) = R(:) - a(nu,:).^2;
        % determine next pivot
        [~,ik] = max(R);
    end
    
end
end

function [s] = returnsim(u,kj,M)
    % scaling factor, should not be choosen too small 
    sigma2 = 0.04;
    
    d = zeros(M,1);
    
    for j = 1:M
        d(j) = norm(u(kj,:)-u(j,:))^2 + eps;
    end

    % avoid zeros!
    s = exp(-d./sigma2);
end


