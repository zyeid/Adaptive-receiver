function [I1] =  attack( I0,n )
%ATTACK [attacked_image,attack_name]=attack(I

%I1=uint8(I0);
switch n
    case 0
        %pas d'attaque
        I1=I0;
        %a_name='sans attaque';
    case 1 
        % 3x3 median filter       
        I1=medfilt2(I0,[3 3]);
        %a_name='3x3 median filter';
    case 2
        % 5x5 median filter        
        I1=medfilt2(I0,[5 5]);
        %a_name='5x5 median filter';      
    case 3
        % color reduction : 8 bits --> 4 bits       
        I1=round(I0)/16;
        %a_name='CR:8->4';
    case 4
        % IS(1/8) : intensity scaling
        I1=I0/8;
        %a_name='IS(1/8)';
    case 5
        % Gamma Correction (1.2)        
        I1= uint8(255.*(double(I0)./255).^1.2);
        %a_name='Gamma 1.2';
    case 6
        I1= uint8(255.*(double(I0)./255).^0.8);
        %a_name='Gamma 0.8';
    case 7
        %JPEG QP=40;
        imwrite(I0,'I_t.jpg','Quality',40);
        I1=imread('I_t.jpg');
       % a_name='JPEG 40';
    case 8
        %JPEG QP=30;
        imwrite(I0,'I_t.jpg','Quality',30);  
        I1=imread('I_t.jpg');
        %a_name='JPEG 30';
    case 9
        %JPEG QP=5;
        imwrite(I0,'I_t.jpg','Quality',20);  
        I1=imread('I_t.jpg');
        %a_name='JPEG 20';        
    case 10
        %poisson 
        I1 = imnoise(I0,'poisson');
        %a_name='poisson';
    case 11
        %gaussian (0/0.05) 
        M =0;
        V=0.05;
        I1 = imnoise(I0,'gaussian',M,V);
        %a_name='gaussian (0/0.05)';   
    case 12
        %salt&peper 0.01
        I1 = imnoise(I0,'salt & pepper',0.01);
        %a_name='salt&peper 0.01';
    case 13
        %histeq
        I1 = histeq(I0);
        %a_name='histeq';
    case 14
        % image sharpening
        H = fspecial('unsharp');
        I1 = imfilter(I0,H,'replicate');
        %a_name='sharpening'; 
    case 15
        h = ones(3,3) / 9;
        I1 = imfilter(I0,h);
        %a_name='average filter 5x5'; 
    case 16
        h = ones(5,5) / 25;
        I1 = imfilter(I0,h);
        %a_name='average filter 7x7';              
    case 17
        I1 = imnoise(I0,'speckle',0.01);
        %a_name='speckle noise 0.01';  

end
I1=uint8(I1);
end


