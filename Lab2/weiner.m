clear, close all;

global PATH
PATH = 'G:/DIP/Lab2/';

I = imread(strcat(PATH,'coins_blurred.tif'));
[row, col] = size(I);

h = fspecial('disk',2);
H = fftshift(fft2(h,row,col));
figure;
mesh(abs(H));
colormap;
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

J = inv_filter(I, H);
figure;
montage([I J], 'Size', [1 1]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

J = wiener_filter(I, H);
figure;
montage([I J], 'Size', [1 1]);
set(gcf, 'units','normalized','outerposition',[0 0 1 1]);



function J = wiener_filter(I, H)
    
    G = fftshift(fft2(I));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    NSR = 0.005;
    F = (conj(H) ./ ((abs(H).^2) + NSR)).*G;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    figure;
    mesh(abs(F));
    colormap;
    title('restored dist. - wnr filter');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    %I2 = (ifft2(ifftshift(F), 'symmetric'));
    I2 = (ifft2(ifftshift(F)));
  
    J = real(I2);
    J = (J <=255).*J + (J > 255).*255;
    J = uint8(J)
end


function J = inv_filter(I, H)
    
    G = fftshift(fft2(I));

    [r, c] = size(G);
    F = zeros(r,c);

    for u = 1:r
        for v = 1:c
            if abs(H(u,v)) > 0.00
                F(u,v) = G(u,v) / H(u,v);
            else
                %F(u,v) = G(u,v);
            end
        end
    end

    figure;
    mesh(abs(F));
    colormap;
    title('restored dist. - inv filter');
     set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    I2 = (ifft2(ifftshift(F)));
  
    J = real(I2);
    J = uint8(J)
    
end