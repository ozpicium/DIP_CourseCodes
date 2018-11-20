clear , close all;

global PATH
PATH = 'G:/DIP/Lab1/';

global IMG
IMG = ["Lena.bmp", "Mandrill.bmp", "Peppers.bmp", "pout.tif", "flower.tif"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I1 = imread(strcat(PATH,char(IMG(5))));
I2 = imread(strcat(PATH,char(IMG(4))));
%I2 = [] %%testing hist. equalization using matching/transfer

histTransfer(I1, I2);

waitforbuttonpress;
close all;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = histTransfer(I1, I2)
    global PATH 
    global IMG
    
    Iorg = I1;
    hst = imhist(Iorg);
    cdf = cumsum(hst) / numel(Iorg); 
    figure;
    subplot(1,3,1);
    imshow(Iorg);
    title('original image');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
  
    if isempty(I2)  %%test for equalization
        hst_ref = ones(1, 256, 'double') * (double(numel(Iorg))/256);
        cdf_ref = cumsum(hst_ref) / numel(Iorg);
    else
        Iref = I2;

        hst_ref = imhist(Iref);
        cdf_ref = cumsum(hst_ref) / numel(Iref);
        
        waitforbuttonpress;
        subplot(1,3,2);
        imshow(Iref);
        title('match image');
        set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    end

    TF = zeros(1,256,'uint8'); %%tranfer function
    
    for i = 1 : 256
        [h,idx_ref] = min(abs(cdf_ref - cdf(i)));
        TF(i) = idx_ref;   %% Tranfer mapping from original image level to reference level
    end

    I2 = TF(Iorg);
    hst_fin = imhist(I2);
    cdf_fin = cumsum(hst_fin) / numel(I2);
    waitforbuttonpress;
    subplot(1,3,3);
    imshow(I2);
    title('result');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);

    waitforbuttonpress;
    x = 1:256;
    figure;
    subplot(1,3,1);
    bar(x, hst);
    xlabel('original hist');
    subplot(1,3,2);
    bar(x, hst_ref);
    xlabel('ref hist');
    subplot(1,3,3);
    bar(x, hst_fin);
    xlabel('final hist');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
    
    waitforbuttonpress;
    figure;
    subplot(1,3,1);
    plot(x, cdf);
    xlabel('original cdf');
    subplot(1,3,2);
    plot(x, cdf_ref);
    xlabel('ref cdf');
    subplot(1,3,3);
    plot(x, cdf_fin);
    xlabel('final cdf');
    set(gcf, 'units','normalized','outerposition',[0 0 1 1]);
end




