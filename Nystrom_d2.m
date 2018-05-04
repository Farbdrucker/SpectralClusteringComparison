function [EVEC,EVAL] = Nystrom_d2(Data,L)
%% Use Nystroem extension for full graph clustering
%% derived from Bertozzi and Flenner 
%%
%% DIFFUSE INTERFACE MODELS ON GRAPHS FOR
%% CLASSIFICATION OF HIGH DIMENSIONAL DATA
%%
%% by Lukas Sanner Summer 2017
%%  INPUT:   Data = (X,Y) coordinates of points in 2d 
%%          L is number of eigenvalues to compute (size of subsystem W_XX)
%%
%%  Output: V       is eigenvector
%%          EVAL    is eigenvalue

% N is Number of data points in 2d
[N,~] = size(Data);

% set a seed to compare results
seed = 1;
rs = RandStream('mt19937ar','Seed',seed);

% random perumation for the sets of X and Y
% Z = X \cup Y and X \cab Y = \varnothing
Z = randperm(rs,N);
X =Z(1:L);
Y =Z(L+1:end);  

%% calculate W_XX and W_XY
% Distance between every dot from X and Y
D_XX = Ny_dist2(Data,X,X);

D_XY = Ny_dist2(Data,X,Y);


%% scale
S = 0.04;

% similarity function exp(-(dist_X - dist_Y)^2/sigma)
W_XX = exp(-D_XX./S) + eps; % add eps to avoid true zeros
W_XY = exp(-D_XY./S) + eps;

W_YX = W_XY';

 
%% normalization of W L_sym
% calculate d_X and d_Y
oneL = ones(L,1);
oneNL = ones(N-L,1);

% degree dX and dY are the sum of the weights
dX = W_XX*oneL + W_XY*oneNL;

% solve lin system for inverse of W_XX
% Use Pinv - Moore Penrose inverse!
% W_XX might be singular
dY = pinv(W_XX)*(W_XY*oneNL);

dY = W_YX*oneL + W_YX * dY;

sX = sqrt(dX);
sY = sqrt(dY);
W_XX = W_XX./(sX*sX');
W_XY = W_XY./(sX*sY');
W_YX = W_XY';
%% calculate orhtonomalized matrix Phi for eigenvectors
[B_X,Gamma,~] = svd(W_XX);
B_Xt = B_X';

S = B_X* Gamma^(-1/2) * B_Xt;

Q = W_XX + S'*(W_XY * W_YX)*S;

[A,Theta,~] = svd(Q);

Phi = [ B_X *Gamma^(1/2); W_YX*B_X*Gamma^(-1/2)]* B_Xt*(A*Theta^(-1/2));


%% Output
% column of Phi are EVECS and 1-Theta(i,i) are [EVAL of L_sym]
EVAL = 1-diag(Theta);
V    = Phi;

II = find(EVAL<0);
EVAL(II) = [];
V(:,II) = [];

[m,n] = size(V);
EVEC = zeros(m,n);
EVEC(Z(1:m),:) = V(1:m,:);
end

function D = Ny_dist2(S,X,Y)
    %% Distance function for dots in R^2
    L = length(X);
    K = length(Y);

    D = zeros(L,K);

    for i=1:L
        for j = 1:K
            D(i,j) = norm(S(X(i),:)-S(Y(j),:))^2+eps;  
        end
    end
end


