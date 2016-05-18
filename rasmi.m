clc; clear all; close all;

delta=158;%158
Nt=512;%156       //    512
N=1024;
DD=2;
MD=5;
%%                  --------    Insertion   ---------------
% read image
I=imread('lena.bmp');
%I=rgb2gray(I);

% direct wavelet transform
[CA1,CH1,CV1,CD1] = dwt2(I,'sym2');
[CA2,CH2,CV2,CD2] = dwt2(CA1,'sym2');
[CA3,CH3,CV3,CD3] = dwt2(CA2,'sym2');
[CA4,CH4,CV4,CD4] = dwt2(CA3,'sym2'); 
%select subband(LL4)
    
     
%generate message
[ messages_mtx,m,M,pos ] = gen_msg( Nt,N );

%perform insertion
CD4(2:33,2:33)=ins_Q(CD4(2:33,2:33),messages_mtx,delta);

%inverse wavelet transform
CA3=idwt2(CA4,CH4,CV4,CD4,'sym2');
CA2=idwt2(CA3,CH3,CV3,CD3,'sym2');
CA1=idwt2(CA2,CH2,CV2,CD2,'sym2');
I_1=idwt2(CA1(1:257,1:257),CH1,CV1,CD1,'sym2');
I_1=uint8(I_1);
% imshow(I_1)
 disp(['PSNR: ' num2str(psnr(I_1,I)) 'dB']) 


 
 
 I_w_a=attack(I_1,4);
 
%%              --------    Detection   ---------------


[CA1_w_a,CH1_w_a,CV1_w_a,CD1_w_a] = dwt2(I_w_a,'sym2');%direct wavelet transform
[CA2_w_a,CH2_w_a,CV2_w_a,CD2_w_a] = dwt2(CA1_w_a,'sym2');
[CA3_w_a,CH3_w_a,CV3_w_a,CD3_w_a] = dwt2(CA2_w_a,'sym2');
[CA4_w_a,CH4_w_a,CV4_w_a,CD4_w_a] = dwt2(CA3_w_a,'sym2'); 
 
 
[ m_host,M_host ] = dec_msg( CD4_w_a(2:33,2:33),pos,Nt,N);%decode message


%[ ratio_M_prog_sub, ratio_M_prog] = progressive( m_host,m,M_host,M );%progressive regressive
mli( m_host,m,M_host,M,1)



