clear , close all;

global PATH
PATH = 'G:/DIP/Lab1/';

global IMG
IMG = ["Lena.bmp", "Mandrill.bmp", "Peppers.bmp", "dome.jpg", "uni.jpg"];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 0
    I = imread(strcat(PATH,char(IMG(1))));
    %figure; imshow(I);

    J = histeq(I);
    figure; imshow(J);
    
    x = 1:256;
    figure;
    subplot(1,2,1);
    bar(x, imhist(I));
    xlabel('original hist');
    subplot(1,2,2);
    bar(x, imhist(J));
    xlabel('final hist');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% LOCAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runLocalHistEq(1, 'Lenna');
% runLocalHistEq(2, 'Mandy');
% runLocalHistEq(3, 'Peps');
% runLocalHistEq(4, 'Dome');
runLocalHistEq(5, 'Uni');


%%%%%%%%%%%%%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function runLocalHistEq(idx, tag)
    global PATH
    global IMG
    
    I = imread(strcat(PATH,char(IMG(idx))));    
    Ipad = refelct_padding(I);
    imwrite(Ipad, strcat(PATH, 'results/', 'reflect', tag, '.bmp'), 'bmp');

    kernel = [21 21];
    J = localhisteq(I, Ipad, kernel);
    imwrite(J, strcat(PATH, 'results/', 'localresult', tag, '21.bmp'), 'bmp');

    kernel = [63 63];
    J = localhisteq(I, Ipad, kernel);
    imwrite(J, strcat(PATH, 'results/', 'localresult', tag, '63.bmp'), 'bmp');

    kernel = [105 105];
    J = localhisteq(I, Ipad, kernel);
    imwrite(J, strcat(PATH, 'results/', 'localresult', tag, '105.bmp'), 'bmp');

    kernel = [189 189];
    J = localhisteq(I, Ipad, kernel);
    imwrite(J, strcat(PATH, 'results/', 'localresult', tag, '189.bmp'), 'bmp');
    
    kernel = [441 441];
    J = localhisteq(I, Ipad, kernel);
    imwrite(J, strcat(PATH, 'results/', 'localresult', tag, '441.bmp'), 'bmp');
end


function Ipad = refelct_padding(I)
    Icorner = imrotate(I, 180);
    Ihor = flip(I, 2);
    Iver = flip(I, 1);
    
    IJ1 = cat(2, Icorner, Iver);
    IJ1 = cat(2, IJ1, Icorner);
    IJ2 = cat(2, Ihor, I);
    IJ2 = cat(2, IJ2, Ihor);
    
    Ijoin = cat(1, IJ1, IJ2);
    Ijoin = cat(1, Ijoin, IJ1);
    
    Ipad = Ijoin;
    
end

function J = localhisteq(I, Ipad, kernel)
    global PATH

    %I = imresize(I, 1/2);
    
    [row, col] = size(I);

    I2 = zeros(row, col, 'uint8');
  
    krow = kernel(1);
    kcol = kernel(2);
            
    pr = 'applying local eq. kernel'
    krow

    for r = 1:row
        for c = 1:col
            local_hist = zeros(1, 256);
           
            %%% apply kernel - pixel at center %%%
            for kr = ceil(-krow/2):floor(krow/2)            
                for kc = ceil(-kcol/2):floor(kcol/2)
                    lev = Ipad(row + (r + kr), col + (c + kc));    %%%collect gray level from neighbourhood                  
                    local_hist(lev+1) =  local_hist(lev+1) + 1;  %%increment bin                   
                    %pause
                end
            end
            %%%%%%%%%%%% kernel done %%%%%%%%%%%%%%%%%%%%%%%
            
            local_cdf = double(cumsum(local_hist))/ (krow*kcol);         
            I2(r,c) = local_cdf(I(r,c)+1) * 255;
        end
    end
    
    J = I2;
end

