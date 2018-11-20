clear, close all;

global PATH
PATH = 'G:/DIP/Lab2/';

global IMG
IMG = ["lena.tif", "flower.tif", "pepper_corrupt.tif", "coins_blurred.tif", "uni.jpg", "dome.jpg", "Mandrill.bmp"];
s = size(IMG);
NI = s(2);

I = imread(strcat(PATH,char(IMG(1))));
figure;
imshow(I);

[row, col] = size(I);
J = ButterWorth(I, 100, 5.0, 1);
waitforbuttonpress;
figure;
imshow(J);

%close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = ButterWorth(I, D0, n, type)  %% 0: low pass;  1: high pass
    sft = 1;   %% center shift
    
    F = fft2(I);
    if sft
        F = fftshift(F);
    end
    [row, col] = size(I);
    
    H = zeros(row, col, 'double');
    I2 = zeros(row, col);
    rs = ceil(row/2)*sft;
    cs = ceil(col/2)*sft;
    
    for u=1:row
        for v=1:col
            if type == 0
                H(u,v) = 1 / ( 1 + (((u-rs)^2 + (v-cs)^2)/(D0^2))^n );
            else
                H(u,v) = 1 / ( 1 + ((D0^2)/((u-rs)^2 + (v-cs)^2))^n );
            end
            
            I2(u,v) = F(u,v) * H(u,v);
        end
    end
    
    waitforbuttonpress;
    figure;
    mesh(H);
    colormap;
    
    waitforbuttonpress;
    figure;
    subplot(1,3,1);
    mesh(abs(F));
    colormap;
    subplot(1,3,2);
    mesh(H);
    colormap;
    subplot(1,3,3);
    mesh(abs(I2));
    colormap;
    
    if sft
        I2 = ifftshift(I2);
    end
    J = uint8(real(ifft2(I2)))
    
end