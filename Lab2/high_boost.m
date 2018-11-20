clear, close all;

global PATH
PATH = 'G:/DIP/Lab2/';

global IMG
IMG = ["lena.tif", "flower.tif", "pepper_corrupt.tif", "coins_blurred.tif", "uni.jpg", "dome.jpg", "Mandrill.bmp"];
s = size(IMG);
NI = s(2);

I = imread(strcat(PATH,char(IMG(7))));
J = HBF(I, 0.7);
figure;
montage([I, J]);





function J = HBF(I, A)
    LF = fspecial('laplacian', 0.1)
    Ilf = imfilter(I, LF);
    
    J = ((A-1) * I) + Ilf;
end