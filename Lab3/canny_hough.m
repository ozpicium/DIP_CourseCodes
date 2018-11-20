clear , close all;

global PATH
PATH = 'G:/DIP/Lab3/';

global IG_THETA

part = 2;

if part == 0
        I = rgb2gray(imread(strcat(PATH,'sniper.jpg')));
        I = edge(I, 'canny', [], 2);
        imshow(I);
        imwrite(I, strcat(PATH, 'sniper_canny2.jpg'));

elseif part == 1
    
        %I = imread(strcat(PATH,'lena.tif'));
        I = rgb2gray(imread(strcat(PATH,'sniper.jpg')));
        %imshow(I);

        for GAUSSIAN_EDGE_FILTERING = 1
        % (1)  %%%%% GAUSSIAN EDGE FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%
                [Ig, Ig_theta] = gauss_filter(I);
                Ig = Ig * (255 / max(Ig(:)));   %%Ig values can extend beyond 255...for purpose of image representation scaling this to 255
                Ig = uint8(Ig);
                figure;
                imshow(Ig);
                
                IG_THETA = Ig_theta;
                waitforbuttonpress;
        end

        
        for NON_MAX_SUPPRESS = 1
        % (2)  %%%%%%%% NON_MAX SUPPRESSION %%%%%%%%%%%%%%%%%%%%%%%%%
        
                relax = 1; %%%%%%  THIS IS THE RELAXATION VALUE USED TO INCLUDE NON_MAX EDGES WHICH ARE STILL CLOSE TO MAX...GIVES SMOOTHNESS
                Inmx = non_max_sup(Ig, Ig_theta, relax);
                Inmx = Inmx * (255 / max(Inmx(:)));     %% Scale Edges
                figure;
                imshow(Inmx, []);

                waitforbuttonpress;
        end


        for HYSTER_THRESH = 1
        % (3)  %%%%%%%% HYSTERISIS THRESHOLDING %%%%%%%%%%%%%%%%%%%%%%%%%

                %%% plotting hist to identify frequency of edges of different intensities....HIGH and LOW threshold are set on the basis of distribution %%%
                ah = imhist(Inmx);
                figure;
                plot(ah(10:255));
                set(gca,'XTick',(0:10:250))

                IH = hys_thresh(Inmx, Ig_theta, 20, 30);
                figure;
                imshow(IH);
                
                waitforbuttonpress;
        end
        
        imwrite(IH, strcat(PATH, 'sniper_canny.jpg'));
        close all;
        %%%%%%%%%%%%%%%%%%%%  PART 1 DONE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


elseif part == 2

        %I = rgb2gray(imread(strcat(PATH,'sniper.jpg')));
        I = imread(strcat(PATH,'sniper_canny.jpg'));
        figure;
        imshow(I);
        
        %%%%%%%%% HOUGH TRNSFORM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [rows, cols] = size(I);
% 
%         [extJ, rho, tht] = hough(I);
%         imshow(I);
%         imshow(imadjust(rescale(H)),'XData',T,'YData',R,...
%             'InitialMagnification','fit');
      
    for MAX_RHO_THETA_EXTRACTION = 1
        if ~MAX_RHO_THETA_EXTRACTION
            break;
         end
        
                
        [accum, rhoSz, thtSz] = huff(I);   %%% get HOUGH %%%%%%%%%%%%%%

        accum = (accum > 0).*accum;           %%%threshold check
        %numel(accum ~=0)
        
        J = uint8(accum);
        J = imresize(J, [rows, cols]);
        figure;
        imshow(J);
            
        %%%%%%%%%%%%% MAX extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [mx, idx] = sort(accum(:), 'descend');
        [rhO thetA] = ind2sub(size(accum), idx(1:1));
        rhO = (rhO - ceil(size(accum, 1)/2)) * rhoSz   %%% until now, only quantization was done...now scaling back to original scale with quant
        thetA = thetA * thtSz - 90                      %%% actual angles are -90 to 90
        
%         rhO = rhO ;%+ 25*rhoSz;%+ 25;
%         thetA = thetA - (thtSz/2);% -2;

    end
        
    for LINE_CALCULATION = 1    
         if ~LINE_CALCULATION
            break;
         end
        
                         
         dispersion = Inf;
         for z = 1:length(rhO) 
                %%%%%  LINE CALCULATION %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                rho = rhO(z);
                theta = thetA(z);
                
                theta = theta * (pi/180);   %%% radians
                %%% equation: rho = x*cos(theta) + y*sin(theta)
                if sin(theta) == 0
                        Ytmp = 1:rows;
                        Xtmp = round(ones(1,rows) * rho / cos(theta));
                elseif cos(theta) == 0
                        Xtmp = 1:cols;
                        Ytmp = round(ones(1,cols) * rho / sin(theta));
                else
                        m = -cos(theta) / sin(theta);          %%slope
                        c = rho / sin(theta);                %% intercept

                        xmin = ceil(-cols/2);
                        xmax = floor(cols/2);
                        X = xmin:xmax;
                        Y = round(m*X + c);
                           
                        Xtmp = X + abs(xmin) + 1;  %%convert back to actual indices
                        Ytmp = Y + abs(ceil(-rows/2)) + 1;                     
                end
                        %%% line making ends
                        %%% the Ytmp values obtained above are not bound to the image
                        %%% dimensions becuase they are in continuous Cartesian plane.
                        %%% Hence, it is necessary to first prune the values which have
                        %%% valud membership to the image. The try-catch method used
                        %%% here is a quick but (kind-of) dirty trick to achieve this.
                        XX = [];
                        YY = [];
                        for yt = 1:length(Ytmp)
                            try
                                I(Ytmp(yt), 1);         %% if Ytmp value is not a valid column index of the image, this will cause 'skip'
                                if I(Ytmp(yt), Xtmp(yt)) ~= 0       %%% it also belongs to an edge
                                    YY = [YY, Ytmp(yt)];
                                    XX = [XX, Xtmp(yt)];
                                end
                            catch
                                %yt = yt+1;
                            end
                        end
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRUNING THE SPAN OF LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                XXt = []; YYt = [];
                xm = mean(XX);
                ym = mean(YY);
                for k = 1:length(XX)
                    euD = sqrt( (XX(k)-xm)^2 + (YY(k) - ym)^2 );
                    if euD <= 400
                        XXt = [XXt, XX(k)];
                        YYt = [YYt, YY(k)];
                    end
                end
                XX = XXt;
                YY = YYt;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                      DRAW LINE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                if ~isempty(XX) && ~isempty(YY) && (length(XX) > 5)
                       %%%%%%%% TRIAL 1 : measure of VARIANCE %%%%%%%%%%%%%%%%%
                       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        XXD = ((XX - mean(XX)).^2) ./ cols;
                        YYD = ((YY - mean(YY)).^2) ./ rows;
                        XYD = (XXD + YYD) .* (sign(XX).*sign(YY)); 
                        
                        disp = var(XYD);           %%variance of points distribution on the line
                          
                        if disp < dispersion
                            dispersion = disp
                            h = imline(gca, [XX(1) XX(length(XX))], [YY(1) YY(length(YY))]);        %% make line only if dispersion is low
                        end
                end

         end  %%% line making ends...really!!

                   
            
            %%%%%%%% BURN LINE %%%%%%%%%%%%%%%%
            I2 = imread(strcat(PATH,'sniper.jpg'));
            
            binImg = h.createMask();
            I2(binImg) = 255;
            figure;
            imshow(I2);
            
    end

    
    for BLOCKWISE_LINE_BUILDING = 0
        if ~BLOCKWISE_LINE_BUILDING
            break;
        end
        
        div = 30;   %%%%%%%%%%%%  <--------------- PUT BLOCK VALUE HERE !!!!!!!!!!!!!!!!!!!!!!!
        
        org = I;
        fin = [];
        
        
        rB = unique(floor((1:rows)./div))*div       %% contains quantized indices starting from 0
        cB = unique(floor((1:cols)./div))*div
        
        for i = 1:length(rB)-1
            fin_stack = [];
           
            for j = 1:length(cB)-1              
                
                I = imcrop(org, [cB(j)+1  rB(i)+1  cB(j+1)-cB(j)  rB(i+1)-rB(i)]);  %%block image
                
                %%%%%%%%% HOUGH TRNSFORM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [rows, cols] = size(I);

                [accum, rhoSz, thtSz] = huff(I, [2 1]);   %%% get HOUGH %%%%%%%%%%%%%%

                accum = (accum > 254).*accum;           %%%threshold check
                %numel(accum ~=0)


                for MAX_RHO_THETA_EXTRACTION = 1
                %%%%%%%%%%%%% MAX extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [mx, idx] = max(accum(:));
                [rho theta] = ind2sub(size(accum), idx);
                rho = (rho - ceil(size(accum, 1)/2)) * rhoSz;   %%% until now, only quantization was done...now scaling back to original scale with quant
                theta = theta * thtSz;

%                 J = uint8(accum);
%                 J = imresize(J, [rows, cols]);
%                 figure;
%                 imshow(J);
                end

                for LINE_CALCULATION = 1        
                %%%%%  LINE CALCULATION %%%%%%%%%%%%%%%%%%%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                theta = theta * (pi/180);   %%% radians
                %%% equation: rho = x*cos(theta) + y*sin(theta)
                if sin(theta) ~= 0
                    m = -cos(theta) / sin(theta);          %%slope
                    c = rho / sin(theta);                %% intercept

                    xmin = ceil(-cols/2);
                    xmax = floor(cols/2);
                    X = xmin:xmax;
                    Y = round(m*X + c);

                    Xtmp = X + abs(xmin) + 1;  %%convert back to positive indices
                    Ytmp = Y + abs(ceil(-rows/2)) + 1;

                    %%% the Ytmp values obtained above are not bound to the image
                    %%% dimensions becuase they are in continuous Cartesian plane.
                    %%% Hence, it is necessary to first prune the values which have
                    %%% valud membership to the image. The try-catch method used
                    %%% here is a quick but (kind-of) dirty trick to achieve this.
                    XX = [];
                    YY = [];
                    for yt = 1:length(Ytmp)
                        try
                            I(Ytmp(yt), 1);         %% if Ytmp value is not a valid column index of the image, this will cause 'skip'
                            YY = [YY, Ytmp(yt)];
                            XX = [XX, Xtmp(yt)];
                            k = k+1;
                        catch
                            %yt = yt+1;
                        end
                    end
                    YY;
                    XX;

                end
                        
                     figure(1);
                     imshow(I);
                     if ~isempty(XX) && ~isempty(YY)
                         %I = ones(size(I,1), size(I,2), 3, 'uint8').*I;
                         h = imline(gca, [XX(1) XX(length(XX))], [YY(1) YY(length(YY))]);
                         binImg = h.createMask();
                         I(binImg) = 255;
                     end
                
                end
                
                if j == 1 
                    fin_stack = I;
                else
                    fin_stack = cat(2, fin_stack, I); 
                end
                
            end
            
            if i == 1 
                    fin = fin_stack;
            else
                    fin = cat(1, fin, fin_stack);  
            end
            figure(2);
            imshow(fin);
            
        end
        
        fin;
        figure;
        subplot(1,2,1);
        imshow(org);
        subplot(1,2,2);
        imshow(fin);
        
    end

end








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%    FUNCTIONS
%%%%%%%%%%%%%%%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Gaussian Edge Filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [J, Jth] = gauss_filter(I, pars)

    if nargin == 2
        sig = pars;
    else
    sig = input('Please enter the sigma value for Gaussian Filter: ');
    end
    
    [row, col] = size(I);
    
    %%%%%%%% KERNEL SIZE CALCULATION %%%%%%%%%%%%%%%%%%%%%
    krow = (2*ceil(sig) + 1);
    kcol = krow;
    
    Ipad = reflect_padding(I);
    
    I2 = zeros(row, col, 'uint8'); %%holder
    
    %%% X, Y edge filtering %%%%
    IX = gausXY(Ipad, sig, krow, kcol, row, col, 'X');
    IY = gausXY(Ipad, sig, krow, kcol, row, col, 'Y');
    
    %%%% Combination %%%%%%%%%%%%%%%
    Ig = sqrt(IX.*IX + IY.*IY);
    
    Ig_theta = round(atan(IY ./ IX) * (180 / pi) / 45) * 45;   %%%%%% conversion to degrees is not necessary...done for convenient visualization!!
    for i = 1:row
        for j = 1:col
            if isnan(Ig_theta(i,j))
                Ig_theta(i,j) = 0;
            end
        end
    end
    
    J = Ig;
    Jth = Ig_theta;
end

%%%%%%% gaussian X & Y filter %%%%%%%%%%%%%%%%%
function J = gausXY(Ipad, sigma, krow, kcol, row, col, dim)
    
%     krow = 5;
%     kcol = 5;
    
    J = zeros(row, col); %%holder
    
    g = zeros(krow,kcol); 
     
    i0 = ceil(krow/2);
    j0 = ceil(kcol/2);
    for i = ceil(-krow/2):floor(krow/2) 
        for j = ceil(-kcol/2):floor(kcol/2)
            ii=i+i0;
            jj=j+j0;
            if dim == 'X'
                g(ii,jj)=(-j / (sigma^2)) * exp( -((i)^2+ (j)^2)/(2*sigma*sigma) );
            else 
                g(ii,jj)=(-i / (sigma^2)) * exp( -((i)^2+ (j)^2)/(2*sigma*sigma) );
            end
        end
    end
    %normalize gaussian filter
    %g = rot90(g,2*90);
    g = -1*g;
    sumg = max(g(:));
    g = g / sumg
    
    
    %%%%% CONVOLUTION %%%%%%%%%%%%%%%%%%
    for r = 1:row
        for c = 1:col
            gaus = 0;
            
            %%% apply kernel - pixel at center %%%
            for kr = ceil(-krow/2):floor(krow/2)            
                for kc = ceil(-kcol/2):floor(kcol/2)
                    lev = double(Ipad(row + (r + kr), col + (c + kc)));                  
                    gaus = gaus + (lev * g(i0+kr, j0+kc));
                end
            end
            %%%%%%%%%%%% kernel done %%%%%%%%%%%%%%%%%%%%%%%
            gaus;
            J(r,c) = gaus;

        end
    end    
 
end

%%%% Reflect Padding %%%%%%%%%%%
function Ipad = reflect_padding(I)
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


%%%% NON MAX SUPPRESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function J = non_max_sup(Ig, Ith, relax)
    
    %%Suppression Kernel Size : 3 x 3
    
    [row, col] = size(Ig);
    
    JJ = zeros(row, col, 'uint8');
    
    %%%% NOTE: The borders are avoided here. Zero padding can be done to
    %%%% include borders. However, it seems redundant given the size of the
    %%%% image matrix and that suppression kernel is only 3 x 3
    
    for i = 2:row-1
        for j = 2:col-1
            
            th = Ith(i,j);      %%% get gradient angle
            if th == 0   %%% horizontal
                mx = max([Ig(i, j-1) Ig(i,j) Ig(i, j+1)]);      %%% Take max of 3 cells
                JJ(i,j) = uint8(abs(mx - Ig(i,j))<=relax)*Ig(i,j);
            elseif th == 45
                mx = max([Ig(i-1, j+1) Ig(i,j) Ig(i+1, j-1)]);      %%% Take max of 3 cells
                JJ(i,j) = uint8(abs(mx - Ig(i,j))<=relax)*Ig(i,j);
            elseif th == -45
                mx = max([Ig(i-1, j-1) Ig(i,j) Ig(i+1, j+1)]);      %%% Take max of 3 cells
                JJ(i,j) = uint8(abs(mx - Ig(i,j))<=relax)*Ig(i,j);
            elseif th == 90 || th==-90
                mx = max([Ig(i-1, j) Ig(i,j) Ig(i+1, j)]);      %%% Take max of 3 cells
                JJ(i,j) = uint8(abs(mx - Ig(i,j))<=relax)*Ig(i,j);
            else
                th, i, j
            end
                    
        end
    end
    
    %%%%%% Now we have JJ which must contain max value edge cells....Ig is already uint8, so no need for conversion %%%%%
    
    J = JJ;
end


%%%%%%%%% HYSTERIS THRESHOLDING %%%%%%%%%%%%%%%%%%%%%%%
function J = hys_thresh(Inmx, Ith, thL, thH)
    
    [row, col] = size(Inmx);
    
    for i = 1:row
        for j = 1:col
            
            if Inmx(i,j) >= thH
                Inmx(i,j) = 255;
                
            elseif Inmx(i,j) < thH && Inmx(i,j) >= thL
                    theta = Ith(i,j);      %%% get gradient angle
                    
                    if theta == 90 || theta == -90             %%% vertical
                        Inmx(i,j) = uint8(Inmx(i, j-1)>=thH || Inmx(i, j+1)>=thH) * 255;   %%horizontal check
                        
                    elseif theta == -45
                        Inmx(i,j) = uint8(Inmx(i-1, j+1)>=thH || Inmx(i+1, j-1)>=thH) * 255;  
                        
                    elseif theta == 45
                        Inmx(i,j) = uint8(Inmx(i-1, j-1)>=thH || Inmx(i+1, j+1)>=thH) * 255;   
                        
                    elseif theta == 0    %%horizontal
                        Inmx(i,j) = uint8(Inmx(i-1, j)>=thH || Inmx(i+1, j)>=thH) * 255;   %%vertical check
                        
                    else
                        theta, i, j
                    end
                     
            else
                Inmx(i,j) = 0;
            end
            
        end
    end
    
    J = Inmx;
end


function [accum, rhoSz, thtSz] = huff(I, pars)
    
     if nargin == 2
             rhoSz = pars(1);
             thtSz = pars(2);
     else
         rhoSz = input('Please enter the bin size for RHO: ');
         thtSz = input('Please enter the bin size for THETA: ');
     end

    global IG_THETA
    
    %rhoSz = 10;  thtSz = 3;
    
    [row, col] = size(I);
    r0 = ceil(row/2);
    c0 = ceil(col/2);
    
    
    thetas = -90:90;
    
    rhoMax = ceil(sqrt(r0^2 + c0^2));     %%% distance from center to extreme vertex
    RHO = -rhoMax : rhoMax;                 %%% not -rhoMax:rhoMax because this is used as index of accumulator 
    
    
    %%%%%%%%%% Quantize as per bin sizes %%%%%%%%%%%%%%%
    THETA = ceil((thetas+90) ./ thtSz);
    THETA = ~THETA + THETA;             %%% to put first 0 into bin 1 (convert 0 to 1)
    RHO = ceil(RHO ./ rhoSz);
    
    accum = zeros(size(unique(RHO), 2), size(unique(THETA), 2));      %%counts -- accumlator matrix
    size(accum);
    
    
    for r = 1:row
        rs = r - r0;        %% shited
        for c = 1:col
            cs = c - c0;
            
            if I(r,c) == 255      %%% it's an edge pixel
                C = sqrt(rs^2 + cs^2);
                ph = atan(cs / rs);
                if isnan(ph)
                    ph = 0;
                end
                
                rhos = C * sin((thetas * (pi/180)) + ph);   %%% vector of rhos
                
                %%%%% GETTING ACCUMULATOR INDICES %%%%%%%%%%%%%%%%%%%
                rhos = ceil(rhos) + rhoMax + 1;       %%% converting to accum RHO index
                rhos = ceil(rhos ./ rhoSz);
                thts = ceil((thetas+90) ./ thtSz);       %%% converting to accum THETA index
                thts = ~thts + thts;            %%% convert zero to one...becuase matrix does not have zero indez
                
                
                for z = 1:size(rhos,2)
                    %rhos(z), thts(z)
                    cond = 0;
                    ornt = round((thts(z)*thtSz - 90) / 45)*45;  %%quantize to sectors
                    
                    %%%% check if line orientation is same as edge orientation at this point
                    if ornt == IG_THETA(r,c)  
                        cond = 1;
                    elseif (abs(ornt) == 90) && (abs(ornt)==abs(IG_THETA(r,c)))
                        cond = 1;
                    end
                    
                    if cond == 1
                        accum(rhos(z), thts(z)) = accum(rhos(z), thts(z)) + 1;
                    end
                end
                
                
            end
        end
    end

    
    accum = accum * (255 / max(accum(:)));
    

end
