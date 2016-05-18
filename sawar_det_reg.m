function [ output_args ] = sawar_det_reg( m_host,m,param_mli,titre )
%SAWAR_DET_REG Sawar detection regions
    T(1,:)=m_host;
    T(2,:)=m;
    T=sortrows(T',2)';   
    T1=T(1,T(2,:)==1);
    T0=T(1,T(2,:)==0);
    
    y=min(m_host):0.01:max(m_host);
    %methode de detection mli ou z
    %y_bit  = det_mli( y,param_mli(1),param_mli(2),param_mli(3),param_mli(4));
    y_bit  = det_z( y,param_mli(1),param_mli(2),param_mli(3));
    
    hold on
        plot(T1,1:size(T1,2),'.','color', 'r')
        plot(T0,1:size(T0,2),'.','color', 'b')
        plot(y,y_bit*max(m_host))
    hold off    

    title(titre)
end

