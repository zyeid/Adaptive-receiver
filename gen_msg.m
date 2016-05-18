function [ messages_mtx,m,M,pos ] = gen_msg( Nt,N )

%Nt= size of pilot message
%N= size of complete message
%messages_mtx: message(m,M) à insérer de taille 32x32
%m: pilot message "ligne"
%M: hiden message "ligne"
%pos: position de m et M dans messages_mtx  1:Nt==>m    //   Nt+1:1024==>M
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Nt=1000;%156


indexes = randperm(Nt);
m = zeros(1, Nt);
m(indexes(1:Nt/2)) = 1; % message connu
M=round(rand(1,N-Nt));% hiden data
    
message_c(1:Nt)=m;
message_c(Nt+1:N)=M;
pos=randperm(N);

message_p(pos)=message_c;% permuted message

%convert to matrix
for i=1:sqrt(N)
    x=(i-1)*sqrt(N)+1;
    messages_mtx(i,1:sqrt(N))=message_p(x:x+sqrt(N)-1);
end


end
