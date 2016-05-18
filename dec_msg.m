function [ m_host,M_host ] = dec_msg( CA4,pos,Nt,N )

%Nt=156;
for i=1:sqrt(N)
    x=(i-1)*sqrt(N)+1;
    CA4_v(x:x+sqrt(N)-1)=CA4(i,:);
end

for i=1:N
    message_r(i)=CA4_v(pos(i));
end
m_host=message_r(1:Nt);
M_host=message_r(Nt+1:N);


end

