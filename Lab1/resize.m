clear, close all
PATH = 'G:/DIP/Lab1/';

IMG = ["Lena.bmp", "Mandrill.bmp", "Peppers.bmp"];
s = size(IMG);
NI = s(2);

quants = [2, 4, 8, 16];
s = size(quants);
NQ = s(2);

for i=1:NI
    
    I = imread(strcat(PATH,char(IMG(i))));
    figure;
    imshow(I);
    title('original');

    for l = 1:NQ
        q = double(quants(l));
        I2 = imresize(I, 1/q);
        I2 = imresize(I2, q);
        waitforbuttonpress;
        figure;
        imshow(I2);
        title(strcat('division factor:', string(q)));
    end
    
    waitforbuttonpress;
    close all;
     
end
