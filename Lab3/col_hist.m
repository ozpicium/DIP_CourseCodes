clear , close all;

global PATH
PATH = 'G:/DIP/Lab3/';

I1 = imread(strcat(PATH,'flowers.jpg'));
I2 = imread(strcat(PATH,'sniper.jpg'));

%I = rgb2gray(I1);
I = I1;
color_hist_eq(I);


function color_hist_eq(I)
    
    org = I;  %%% backup
    
%%%%%% RGB eq %%%%%%%%%%%%%%%%%%%
    R = I(:, :, 1);
    G = I(:, :, 2);
    B = I(:, :, 3);

    R = my_hist_eq(R);
    G = my_hist_eq(G);
    B = my_hist_eq(B);
    
    J1 = cat(3, R, G, B);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% HSV eq %%%%%%%%%%%%%%%%%%%
    I = rgb2hsv(org);
    
    H = I(:, :, 1);
    S = I(:, :, 2);
    V = I(:, :, 3)*255;

    V = my_hist_eq(V)/255;

    J2 = cat(3, H, S, V);
    J2 = hsv2rgb(J2)*255;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Lab eq %%%%%%%%%%%%%%%%%%%
    I = rgb2lab(org);
    L = I(:, :, 1);
    a = I(:, :, 2);
    b = I(:, :, 3);

    L = my_hist_eq(L, 101);
    %L = histeq(L)
    
    J3 = cat(3, L, a, b);
    J3 = lab2rgb(J3, 'OutputType','uint8');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    figure;
    montage([org J1 J2 J3], 'Size', [1 1]);
    title('ORG.      RGB.      V.      L.     ');
    
%     J1 = my_hist_eq(G);
%     J2 = histeq(G);
%     figure;
%     montage([G J1 J2], 'Size', [1 1]);
%     
%     J = J1;
end


function J = my_hist_eq(I, n)
    if nargin == 1
        n = 256;
    end
    
    if 0
        hgram = ones(1,n);
        hgram = (numel(I)/n)*hgram;
        J = histeq(I, hgram);
    else
        hist = zeros(1, n);
        [r, c] = size(I);

        J = zeros(r, c);

        for i = 1:r
            for j = 1:c
                try
                    lev = round(I(i,j));
                    hist(lev+1) =  hist(lev+1) + 1;  %%%% increment bin count
                catch
                    lev
                end
            end
        end


    %     cdf = 1:n;
    %     cdf = cdf / n;
        cdf = double(cumsum(hist))/ numel(I);
        nHist = (double(hist)/numel(I));
        nHist = nHist*(0.5/max(nHist));
        figure;
        x = 1:n;
        plot(x, cdf);
        hold on
        bar(x, nHist);
        hold on

        for i = 1:r
            for j = 1:c
                lev = round(I(i,j));
                J(i,j) = n * cdf(lev+1);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% plotting final cdf for verification %%%%%%%%%%%%%%%%%%
        hist = zeros(1, n);
        [r, c] = size(J);

        for i = 1:r
            for j = 1:c
                try
                    lev = round(J(i,j));
                    hist(lev+1) =  hist(lev+1) + 1;  %%%% increment bin count
                catch
                    %lev
                end
            end
        end
        cdf = double(cumsum(hist))/ numel(J);
        plot(x, cdf);
        xlim([0, n]);
        set(gca,'XTick',(0:10:n))
        hold off
        
    end
end