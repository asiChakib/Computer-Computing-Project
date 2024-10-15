%-------------------------------------------------------
% Calcul des vecteurs propres et les valeurs propre 
% en se basant sur subspace_iter_v0 
%-------------------------------------------------------

%%  Application de la SVD : compression d'images

clear all
close all

% Lecture de l'image
I = imread('BD_Asterix_0.png');
I = rgb2gray(I);
I = double(I);
[q, p] = size(I);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% On choisit de ne considérer que 200 vecteurs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QUELQUES VALEURS PAR DÉFAUT DE PARAMÈTRES, 
% VALEURS QUE VOUS POUVEZ/DEVEZ FAIRE ÉVOLUER

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tolérance
eps = 1e-8;

% nombre d'itérations max pour atteindre la convergence
maxit = 10000;

% taille de l'espace de recherche (m)
search_space = 400;

% pourcentage que l'on se fixe
percentage = 0.995;

% p pour les versions 2 et 3 (attention p déjà utilisé comme taille)
puiss = 1;

%%%%%%%%%%%%%
% À COMPLÉTER
%%%%%%%%%%%%%

%%
% calcul des couples propres
%%
if p<q
    I_C = I*I';
else 
    I_C = I'*I;
end 

fprintf('subspace v1\n')
tic
[ W, D, n_ev, it, itv, flag ] = subspace_iter_v1( I_C, search_space, percentage, eps, maxit );
toc 

size(W)
D = diag(D);    


%%
% calcul des valeurs singulières
%%
sigmas = sqrt(sum(D,2));

%%
% calcul de l'autre ensemble de vecteurs
%%

V= zeros(p,p);
U= zeros(q,q);    

if p<q 
    U = W;
    for k =1:size(W,2)
        V(:,k)= (1/sigmas(k))*I'*W(:,k);
    end
else
    V = W;
    for k =1:size(W,2)
        U(:,k)= (1/sigmas(k))*I*W(:,k);
    end
end    





% vecteur pour stocker la différence entre l'image et l'image reconstuite
inter = 1:30:size(W,2);
differenceSVD = zeros(size(inter,2), 1);

%%
% calcul des meilleures approximations de rang faible
%%
ti = 0;
td = 0;
sigmas = diag(sigmas);
for k = inter
    
    % Calcul de l'image de rang k
    Im_k = U(:, 1:k)*sigmas(1:k, 1:k)*V(:, 1:k)';
    % Affichage de l'image reconstruite
    ti = ti+1;
    figure(ti)
    colormap('gray')
    imagesc(Im_k)
    
    % Calcul de la différence entre les 2 images
    td = td + 1;
    differenceSVD(td) = sqrt(sum(sum((I-Im_k).^2)));
    pause
end

% Figure des différences entre image réelle et image reconstruite
ti = ti+1;
figure(ti)
hold on 
plot(inter, differenceSVD, 'rx')
ylabel('RMSE')
xlabel('rank k')
pause