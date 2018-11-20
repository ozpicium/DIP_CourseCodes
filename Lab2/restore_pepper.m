clear, close all;

global PATH
PATH = 'G:/DIP/Lab2/';

I = imread(strcat(PATH,'pepper_corrupt.tif'));

% I2 = medfilt2(I, [3 3]);
% I3 = I - I2;
% I3 = medfilt2(I, [10 10]) + I3;
% montage([I I3], 'Size', [1 1]);
% 
% F = fft2(I);
% F = fftshift(F);
% figure;
% mesh(abs(F));
% colormap;

DL = 20;
DH = 70;
[J,H] = ButterWorth(I, DL, 4, DH, 4, 3);
%J = imsharpen(J, 'radius', 2, 'amount', 8);
%waitforbuttonpress;
figure;
mesh(H);
colormap;

% waitforbuttonpress;
figure;
montage([I J], 'Size', [1 1]);
%imshow(J, [0 max(J(:))]);

O = imread(strcat(PATH,'pepper_corrupt.tif'));
nsr = 0.01;
[Jr, Jc] = size(J);
NSR = zeros(Jr, Jc);

J2 = wiener_filter(J, H, NSR);
%J2 = imsharpen(J2, 'radius', 2, 'amount', 5);
figure;
montage([I J2], 'Size', [1 1]);


function [J,H] = ButterWorth(I, DL, nL, DH, nH, type)  
    sft = 1;   %% center shift
    
    F = fft2(I);
    if sft
        F = fftshift(F);
    end
%     figure;
%     mesh(abs(F));
%     colormap;

    [row, col] = size(I)
    
    H = zeros(row, col, 'double');
    I2 = zeros(row, col);
    rs = ceil(row/2)*sft;
    cs = ceil(col/2)*sft;
    
    for u=1:row
        for v=1:col
            if type == 3
                H(u,v) = (1 / ( 1 + (((u-rs)^2 + (v-cs)^2)/(DL^2))^nL )) + (1 / ( 1 + ((DH^2)/((u-rs)^2 + (v-cs)^2))^nH ));
                I2(u,v) = F(u,v) * H(u,v);
                
            else
                if type == 4
                    H(u,v) = (1 / ( 1 + (((u-rs)^2 + (v-cs)^2)/(DL^2))^nL ));
                    Ft(u,v) = F(u,v) * H(u,v);
                    H(u,v) = (1 / ( 1 + ((DH^2)/((u-rs)^2 + (v-cs)^2))^nH ));
                    I2(u,v) = Ft(u,v) * H(u,v);
                end
            end
        end
    end
    
   
%     
%     waitforbuttonpress;
%     figure;
%     subplot(1,3,1);
%     mesh(abs(F));
%     colormap;
%     subplot(1,3,2);
%     mesh(H);
%     colormap;
%     subplot(1,3,3);
%     mesh(abs(I2));
%     colormap;
    
    if sft
        I2 = ifftshift(I2);
    end
    J = uint8(real(ifft2(I2)));
    
end


function J = wiener_filter(I, H, NSR)
    
    G = fftshift(fft2(I));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = (conj(H) ./ ((abs(H).^2) + NSR)).*G;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     figure;
%     mesh(abs(F));
%     colormap;
%     title('restored dist. - wnr filter');
%     set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    %I2 = (ifft2(ifftshift(F), 'symmetric'));
    I2 = (ifft2(ifftshift(F)));
  
    J = real(I2);
    J = (J <=255).*J + (J > 255).*255;
    J = uint8(J);
end
