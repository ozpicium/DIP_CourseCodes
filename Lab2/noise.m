clear , close all;

global PATH
PATH = 'G:/DIP/Lab2/';

global IMG
IMG = ["lena.tif", "flower.tif"];
s = size(IMG);
NI = s(2);

for i=1:NI
    I = imread(strcat(PATH,char(IMG(i))));
    
    Izmg = imnoise(I, 'gaussian', 0, 0.01);
    imwrite(Izmg, strcat(PATH, 'pic', num2str(i), '_gaussianNoise.png'));
    
    Izmgmf = medfilt2(Izmg, [3 3]);
    MSEmf = immse(Izmgmf, I)
    
    Izmggf = imgaussfilt(Izmg, 1.0);
    MSEgf = immse(Izmggf, I)
    
    figure;
    montage([[I Izmg]; [I Izmgmf]; [I Izmggf]], 'Size', [1 1]);
    
    waitforbuttonpress; %%%%%%%%%%%%%
    
    Isnp = imnoise(I, 'salt & pepper', 0.1);
    imwrite(Isnp, strcat(PATH, 'pic', num2str(i), '_saltnpepperNoise.png'));
     
    Isnpmf = medfilt2(Isnp, [5 5]);
    MSEmf = immse(Isnpmf, I)
    
    Isnpgf = imgaussfilt(Isnp, 1.0);
    MSEgf = immse(Isnpgf, I)
    
    figure;
    montage([[I Isnp]; [I Isnpmf]; [I Isnpgf]], 'Size', [1 1]);
    
    waitforbuttonpress;
end

close all;