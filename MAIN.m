%% TEST DATA
%% aim is to plot:
%%      + first    
%%          - x axis: number of nodes
%%          - y axis: CPU time from EIGS, Nyström, ACA
%%
%%      + second
%%          - x axis: number of eval
%%          - y axis: eigenvalue from EIGS, Nyström, ACA
%%
%% test data derives from 
%%      1 use images with different resolution (to do, due to comparison dots in 2 D are more ordinary to generate)
%%      2 use test samples from crescentfullmoon
%%
%%
%% used functions:
%%      moon_data
%%      matlab_eigs_d2
%%      Nystrom_d2
%%          Ny_dist2
%%          perma
%%      acaf_d2
%%      lanczos_aca
cla 

clear all
close all

%% Data
 total = 3; 
 MinNode = 100;
 NodeTick = 100;
   MaxNode = MinNode+total*NodeTick;
 Nodennz = MinNode + (0:total-1)*NodeTick;


%% Begin testing
%% matlab eigs; Nyström; Adaptive Cross Approximation with Lanczos
% Number of eigenvalues to compute
evalnz = 20;
% for all data sets try functions and quantify the time
functionnz = 3;

TIME = zeros(total,functionnz);
% if any function takes longer then one minute cancel this routine for
% further computation
max_time = 120;

EVAL_EIGS   = zeros(total,evalnz);
EVAL_NY     = zeros(total,evalnz);
EVAL_ACA    = zeros(total,evalnz);

nysubdim    = zeros(total,1);
acarank     = zeros(total,1);

%% Begin loop
for i = 1:total
    fprintf('\n\n%d. test data start\t(%2.1f%%)\n',i,i/total*100)
    fprintf('==================================================\n')
    %% get data
    N = Nodennz(i);
    fprintf('generate data with nodes = %d\n',N)
    X = moon_data(N);

    
    %% eigs
    %----------------------------------------------------------
    j = 1;
    if (i>1 && TIME(i-1,j) > max_time)
        TIME(i,j) = TIME(i-1,j);
        fprintf('eigs was canceld\n')
    end
    fprintf('start eigs')
    tic;
        [Veigs,eval] = matlab_eigs_2d(X,evalnz);
    ttt = toc;
    fprintf('\t--> succeded \t %d s\n',ttt)
    TIME(i,j) = ttt;
    EVAL_EIGS(i,:) =1- diag(eval);

    %% Nyström
    %----------------------------------------------------------
    j = 2;
    % subdimension for nyström method
    %   total dim N and evalnz eigenvalues
    NySubDim = ceil(2*evalnz);
    if NySubDim<2/100*N
        NySubDim = ceil(2/100 * N);
    end
    nysubdim(i) = NySubDim;
    if (i>1 && TIME(i-1,j) > max_time)
        TIME(i,j) = TIME(i-1,j);
        fprintf('nystrom was canceld\n')
    end
    fprintf('start Nystrom')
    tic;
        [Vny,eval] = Nystrom_d2(X,2*evalnz);
    ttt = toc;
    fprintf('\t--> succeded \t %d s\t(%d)\n',ttt,NySubDim)
    Vny = Vny(:,1:evalnz);
    TIME(i,j) = ttt;
    EVAL_NY(i,:) = eval(1:evalnz);

    %% ACA
    %----------------------------------------------------------
    j = 3;
    if (i>1 && TIME(i-1,j) > max_time)
        TIME(i,j) = TIME(i-1,j);
        fprintf('aca was canceld\n')
    end
    fprintf('start ACA')
    tic;
    %[w,acarank(i)] = aca_sym2(X, NySubDim);
    [Vaca, eval] = main_aca(X,evalnz, 2*evalnz);
            
    
    ttt = toc;

   eval =  sort(eval,'descend');
    fprintf('\t--> succeded \t %d s\t(rank %d)\n',ttt,acarank(i))
    TIME(i,j) = ttt;
    for k =1:evalnz
        if k> max(size(eval))
            break;
        end
        EVAL_ACA(i,k) = 1-eval(k);
    end

end

%% cluster
%k = pick_vector(Veigs(:,1:9));
k= 2;
IIeigs = find(Veigs(:,k)<0);
JJeigs = find(Veigs(:,k)>=0);
fprintf('Chosen vector for eigs %d\n',k)

%k = pick_vector(Vaca(:,1:9));
IIaca = find(Vaca(:,k)<0);
JJaca = find(Vaca(:,k)>=0);
fprintf('Chosen vector for ACA %d\n',k)

%k = pick_vector(Vny(:,1:9));
IIny = find(Vny(:,k)<0);
JJny = find(Vny(:,k)>=0);
fprintf('Chosen vector for Nystrom %d\n',k)


figure(1)
subplot(1,3,1)
plot(X(IIeigs,1),X(IIeigs,2),'rx')
hold on
plot(X(JJeigs,1),X(JJeigs,2),'bo')
title('eigs')

subplot(1,3,2)
plot(X(IIny,1),X(IIny,2),'rx')
hold on
plot(X(JJny,1),X(JJny,2),'bo')
title('Nystrom')

subplot(1,3,3)
plot(X(IIaca,1),X(IIaca,2),'rx')
hold on
plot(X(JJaca,1),X(JJaca,2),'bo')
title('ACA')

%% plot data;
%------------------------------------
 figure(3)
 subplot(3,1,1)
 plot(EVAL_EIGS(1,:),'xr')
 hold on
 plot(EVAL_NY(1,:),'ob')
 plot(EVAL_ACA(1,:),'^k')
 xlabel('eigenvalues')
 ylabel('value')
 legend('eigs','Nystrom','ACA','Location','Southeast')
 nameone = sprintf('number of nodes %d',Nodennz(1));
 title(nameone)
 
 subplot(3,1,2)
 mid = ceil(total/2);
 plot(EVAL_EIGS(mid,:),'xr')
 hold on
 plot(EVAL_NY(mid,:),'ob')
 plot(EVAL_ACA(mid,:),'^k')

 %legend('eigs','Nystrom','ACA','Location','Southeast')
 namemid = sprintf('number of nodes %d',Nodennz(mid));
 title(namemid)
 
 subplot(3,1,3)
 plot(EVAL_EIGS(end,:),'xr')
 hold on
 plot(EVAL_NY(end,:),'ob')
 plot(EVAL_ACA(end,:),'^k')
 %legend('eigs','Nystrom','ACA','Location','Southeast')
 nameend = sprintf('number of nodes %d',Nodennz(end));
 title(nameend)
 
 