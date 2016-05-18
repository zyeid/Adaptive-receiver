function [ ratio_M ] = mli( m_host,m,M_host,M,nf_in)

    param_z_pso=PSO_mli(m_host,m);

    
    m_det=det_mli( m_host,param_z_pso(1),param_z_pso(2),param_z_pso(3));
    [number_m,ratio_m] = biterr(m,m_det);
    
    M_det  = det_mli( M_host,param_z_pso(1),param_z_pso(2),param_z_pso(3));
    if ratio_m>0.5
        M_det=not(M_det);
    end
    [number_M,ratio_M] = biterr(M,M_det);
    
   
%     
    % sawar m_host
     T(1,:)=m_host;
     T(2,:)=m;
  
    T=sortrows(T',2)';
   
    T1=T(1,T(2,:)==1);
    T0=T(1,T(2,:)==0);

    
    figure(1)
        title('MLI')
% 
y=min(m_host):0.01:max(m_host);
y_bit  = det_mli( y,param_z_pso(1),param_z_pso(2),param_z_pso(3));
subplot(2,1,1)
  hold on
  
    plot(T1,1:size(T1,2),'.','color', 'r')
    plot(T0,1:size(T0,2),'.','color', 'b')
     plot(y,y_bit*550)%*max(m_host)
    hold off
    title('MLI pilot message')
    % sawar M_host
    T(1,:)=M_host;
     T(2,:)=M;
  
    T=sortrows(T',2)';
   
    T1=T(1,T(2,:)==1);
    T0=T(1,T(2,:)==0);
axis([min(M_host) max(M_host) 0 600])
    
    subplot(2,1,2)

y=min(M_host):0.01:max(M_host);
y_bit  = det_mli( y,param_z_pso(1),param_z_pso(2),param_z_pso(3));
  hold on
    plot(T1,1:size(T1,2),'.','color', 'r')
    plot(T0,1:size(T0,2),'.','color', 'b')
     plot(y,y_bit*550)% *max(M_host)
    hold off
    axis([min(M_host) max(M_host) 0 600])
    title('MLI payload message')
%     
%     
%     
%     if ratio_m>0.5
%         ratio_m=1-ratio_m;
%         number_m=size(m_host,2)-number_m;
%     end
%     disp(['                       DETECTION USING MLI'] )
%     disp(['MLI Pilot message BER: ' num2str(ratio_m) ' , ' num2str(number_m) ' false bit' ] )
%     disp(['MLI Payload message BER: ' num2str(ratio_M) ' , ' num2str(number_M) ' false bit' ] )
end

 function [A]=PSO_mli(LL4,m)
pops=70;%nombre de particule dans la population 70
maxgen=300;% nombre d'it�ration maximal 500

n=3; %n est le nombre des variables de la particule ici deulement delta!
bound=zeros(n,2);%initialisation de la matrice des  limites des variables

% *************Limites des varaiables************
bound(1,:)=[0. 1];%D     % limite de chaque variable  0.15
bound(2,:)=[0 1];%d  0.006
bound(3,:)=[0 2*pi];%gamma   min(min(LL4)) max(max(LL4))


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



end

function [JT]= fitness1(m_host,m,pop)
for i=1:size(pop,1) 
    y=min(m_host):max(m_host);
    if conf_mli(y,pop(i,1),pop(i,2),pop(i,3))==0 
        m_det=det_mli( m_host,pop(i,1),pop(i,2),pop(i,3) );
        [number,ratio] = biterr(m,m_det);
        JT(i)=0.5-abs(ratio-0.5);
    else 
        JT(i)=0.5;
    end
end


end

function[S]=conf_mli(y,F1,F2,C)
y_bit=0.75*5.*sin(2*pi*F2*y+C)>=5.*sawtooth(2*pi*F1*y);
%determiner latgeur de chaque zone
k=1;

for i=1: size(y,2)-1
   if y_bit(i)~=y_bit(i+1) 
       T(k)=y(i);
       k=k+1;
   end
end
if k>=3
k=1;
for i=1:size(T,2)-1
    L(k)=T(i+1)-T(i);
    k=k+1;
end
Mx=max(L)-min(L);
%Me=mean(L);
if Mx>=0.3*max(L)%0.3
    S=1;
else S=0;
end

% if (Mx-Me)>0*Me
%     S=1;
% else S=0;
% end
else
    S=1;
end
end



function [ y_bit ] = det_mli( y,F1,F2,C)

y_bit=0.75*5.*sin(2*pi*F2*y+C)>=5.*sawtooth(2*pi*F1*y);

end




