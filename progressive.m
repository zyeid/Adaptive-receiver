function [ ratio_M2,ratio_M ] = progressive( m_host,m,M_host,M )


 %%figure (1) detection progressive sans subsampling
 DD=2;
 MD=5;
 [ m_host2,m2 ] = reduction( m_host,m,DD,MD );
 param_z_pso2=PSO_z(m_host2,m2);
 m_det2=det_z( m_host2,param_z_pso2(1),param_z_pso2(2),param_z_pso2(3));
 [number_m2,ratio_m2] = biterr(m2,m_det2);
 
  M_det  = det_z( M_host,param_z_pso2(1),param_z_pso2(2),param_z_pso2(3));
  if ratio_m2>0.5
      M_det=not(M_det);
  end
  [number_M2,ratio_M2] = biterr(M,M_det);
    
 figure(1)
 ax1=subplot(3,1,1);%RF
 sawar_det_reg( m_host,m,param_z_pso2,'RF' )
 ax2=subplot(3,1,2);%RF sabsampled
 sawar_det_reg( m_host2,m2,param_z_pso2,'RF subsampled' )
 ax3=subplot(3,1,3);%PL
 sawar_det_reg( M_host,M,param_z_pso2,'PL' )
 linkaxes([ax1,ax2,ax3],'x')
 
 disp(['BER of PL with subsampled method' num2str(ratio_M2)]) 

 
 %%figure(2) detection progressive avec subsampling
 param_z_pso=PSO_z(m_host,m);
 m_det=det_z( m_host,param_z_pso(1),param_z_pso(2),param_z_pso(3));
 [number_m,ratio_m] = biterr(m,m_det);
 
  M_det  = det_z( M_host,param_z_pso(1),param_z_pso(2),param_z_pso(3));
  if ratio_m>0.5
      M_det=not(M_det);
  end
  [number_M,ratio_M] = biterr(M,M_det);
  
 
 figure(2)
 ax1=subplot(2,1,1);%RF
 sawar_det_reg( m_host,m,param_z_pso,'RF' )
 ax2=subplot(2,1,2);%PL
 sawar_det_reg( M_host,M,param_z_pso,'PL' )
 linkaxes([ax1,ax2],'x')
 disp(['BER of PL with Non subsampled method: ' num2str(ratio_M)]) 
 


end

 function [A]=PSO_z(LL4,m)
pops=10;%nombre de particule dans la population 
maxgen=20;% nombre d'it�ration maximal 300

n=3; %n est le nombre des variables de la particule ici deulement delta!
bound=zeros(n,2);%initialisation de la matrice des  limites des variables

% *************Limites des varaiables************
bound(1,:)=[-120 120];%D     % limite de chaque variable 0.5 120
bound(2,:)=[-100 100];%d -2 2
bound(3,:)=[-1000 1000];%gamma   min(min(LL4)) max(max(LL4))


%...
%*************************************************

numvar=n; % Nombre de variables (nombre de dimentions de populations)

rng=(bound(:,2)-bound(:,1))';%largeur de la plage de variation des variables
popi=zeros(pops,numvar);       %popi = initialisation de la population
popi(:,1:numvar)=(ones(pops,1)*rng).*(rand(pops,numvar))+...
    (ones(pops,1)*bound(:,1)');%g�n�ration des valeurs al�atoires de la polulation

wmax=0.9;
wmin=0.4;


c1=1.5;
c2=1.5;


for iter=1:maxgen
W(iter)=wmax-((wmax-wmin)/maxgen)*iter;
end
%on peut �liminer cette boucle et on fixe w � la valeur moyenne

%****************************
popa=popi;
vilocity=rand(numvar,pops)';%initialisation de la vitesse � une valeur al�atoire
Z=popi;
Z2=fitness1(LL4,m,popi);
% for j=2:n
%  for i=1:length(vilocity)
%     if vilocity(i,j-1)>vilocity(i,j);
%    a=vilocity(i,j-1);
%     vilocity(i,j)=a;
%     vilocity(i,j-1)=vilocity(i,j);
%     end
%   end
% end
 
% cette boucle n'est pas n�cessaire elle est utulis� pour un cas pr�cis

k=2; % compteur pour la boucle while

fiti=fitness1(LL4,m,popa); %appel � la fonction � minimiser

[gbest index]=min(fiti);% g�n�ration de la "best postion" global et ses coordonn�es 
xbbest=popa(index,:);% les variables correspondants � la "best position"
ff1(1)=(1/gbest)-1;%d�normalisation de la fonction � minimiser

for i=1:pops
pbest(i)=fiti(i);%g�n�ration de la meilleurs positions de chaque particule de la polpulation
xbest(i,:)=popa(i,:);%g�n�ration des coordonn�es des meilleurs position
end
% verifier l'interet de la boucle %%%%%%%%%%%%%

max_LL4=max(max(LL4));
min_LL4=min(min(LL4));
%%%%%%adaptation%%%%%%
while k<maxgen+1
    i=1:pops;
    vilocity(i,:)=W(k)*vilocity(i,:) + c1*rand(1)*(xbbest(ones(size(popa,1),1),:)-popa(i,:))+c2*rand(1)*(xbest(i,:)-popa(i,:));
    
%% saturation vilocity
for j=1:n
    vilocity(vilocity(:,j) >(bound(j,2)-bound(j,1))/2,j)=(bound(j,2)-bound(j,1))/2;
end

popa(i,:)=popa(i,:)+vilocity(i,:);%  

%% Saturation position
for j=1:n
    popa(popa(:,j) >bound(j,2),j)=bound(j,2);
    popa(popa(:,j) <bound(j,1),j)=bound(j,1);
end
popa(popa(:,2) ==0,2)=0.0000001;


    fitp=fitness1(LL4,m,popa);%appel de la fonction � minimiser
for i=1:pops
    if fitp(i)<pbest(i)
        pbest(i)=fitp(i);%mise � jour de la meilleur position de chaque particile
        
        xbest(i,:)=popa(i,:);%mise � jour des coordonn�es des best
       
    end
end
[gbest index]=min(pbest);%mise � jour de best position globale
xbbest=popa(index,:);% ces coordonn�es xbest

cc(k)=gbest;
ff1(k)=(1/cc(k))-1;% d�normalisation de le fonction � minimiser
k=k+1;
Z=[Z,xbest];
Z2=[Z2,fitness1(LL4,m,xbest)];
end

A=xbbest;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%
% figure(1),
% plot(cc,'b'), hold on;
% xlabel('generation');
% ylabel('fitness');
% title('fitness preogress');
% legend('PSO');

end


function [JT]= fitness1(LL4,m,pop)


for i=1:size(pop,1)

    m_det=det_z( LL4,pop(i,1),pop(i,2),pop(i,3) );
    

%     if isnan(m_det)
%         ratio=0.5;
%     else
    [number,ratio] = biterr(m,m_det);
%     end
    JT(i)=0.5-abs(ratio/4-0.5);
end

%JT
end



function [ wb ] = det_z( lambda,D,d,gama )
%gama=min(min(lambda));
%Detection progressive ou regressive
    S=(-D+sqrt(D^2+2*d.*abs(lambda-gama)))/d;
    wb=mod(floor(abs(S)+0.5),2);%p=0.5
end


