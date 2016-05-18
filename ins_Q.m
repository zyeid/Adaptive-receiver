function [ C1 ] = ins_Q( C,x,S )%     S=delta  x=m
C=double(C);
d(1)=1; d(2)=-1;
C1=q_f(C+sign(d(x+1))*0.25*S,S)-sign(d(x+1))*0.25*S;
end

function [ qC ] = q_f(C, S )
    k=floor(C./S);
    qC=k.*S+0.5*S;
end
