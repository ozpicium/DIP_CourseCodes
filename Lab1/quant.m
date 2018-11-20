clear, close all
PATH = 'G:/DIP/Lab1/';

IMG = ["Lena.bmp", "Mandrill.bmp", "Peppers.bmp"];
s = size(IMG);
NI = s(2);

quants = [2, 4, 8, 16, 64];
s = size(quants);
NQ = s(2);


for i=1:NI
    
    I = imread(strcat(PATH,char(IMG(i))));
    figure;
    imshow(I);
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    title('original');

    for l = 1:NQ
        q = quants(l);
        I2 = double(I) / q;
        I2 = uint8(I2);
        waitforbuttonpress;

        SEL = 1;
        figure;
        if SEL == 1      %%scale to original range
           I2 = I2 * q;
           subplot(1,2,1);
           imshow(I2);
           title(strcat('division factor:', string(q)));
           subplot(1,2,2);
           imshow(I-I2, []);
           set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        else
           %subplot(1,2,1);
           h = imshow(I2, [0 max(I2(:))]);   %%min - max scaling
           title(strcat('division factor:', string(q)));
%            subplot(1,2,2);
%            I2 = I2 * (255 / max(I2(:)));
%            imshow(I-I2, []);
           set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
        end

    end
    
    waitforbuttonpress;
    close all;
     
end
