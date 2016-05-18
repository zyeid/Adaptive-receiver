function [ m_host2,m2 ] = reduction( m_host,m,DD,MD )
%REDUCTION of samples using drop distance(DD%) and merge distance(MD%)
%%              DD reduction
T(1,:)=m_host;
T(2,:)=m; 
TT=sortrows(T',1);

Tx=TT;
i=2;
total=size(TT,1);
dp=DD*(max(m_host)-min(m_host))/100;
while i<total
   if (Tx(i,2)~= Tx(i+1,2)) && (Tx(i+1,1)-Tx(i,1))<dp% même bit
       Tx(i,:)=[];
       Tx(i,:)=[];       
       total=total-2;
       if i>1
           i=i-1;
       end
   else
       i=i+1;
   end    
end
Tx=Tx';
T1=Tx(1,Tx(2,:)==1);
T0=Tx(1,Tx(2,:)==0);

m_host1=[T1,T0];
m1=zeros(1,size(m_host1,2));
m1(1,1:size(m_host1,2)/2)=1;

%%         MD reduction
T_(1,:)=m_host1;
T_(2,:)=m1; 
TT_=sortrows(T_',1);


T_deux=TT_;

i=2;
total=size(T_deux,1);
mp=MD*(max(m_host)-min(m_host))/100;
while i<total
   if(T_deux(i,2)== T_deux(i-1,2)) && (T_deux(i,2)== T_deux(i+1,2)) && (T_deux(i+1,1)-T_deux(i,1))<mp% même bit       
       T_deux(i,:)=[];            
       total=total-1;
       if i>2
           i=i-1;
       end
   else
       i=i+1;
   end    
end

Tx=T_deux';
T1=Tx(1,Tx(2,:)==1);
T0=Tx(1,Tx(2,:)==0);

m_host2=[T1,T0];
m2=zeros(1,size(m_host2,2));
m2(1,1:size(T1,2))=1;



end

