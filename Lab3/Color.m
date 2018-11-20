clear , close all;

global PATH
PATH = 'G:/DIP/Lab3/';

I = imread(strcat(PATH,'flowers.jpg'));
size(I(:))
sz = size(I)
mk = imread(strcat(PATH,'mask.tif'));
mk = uint8(mk);

% sm = 0;
% for i = 1:sz(1)
%     for j = 1: sz(2)
%        if mk(i,j) == 1
%            sm = sm + 1;
%        end
%     end
% end
% sm

%%%%%% SWAP %%%%%%%%%%%%%%%
J1 = swap_RB(I, mk);
figure;
montage([I J1], 'Size', [1 1]);

waitforbuttonpress;

%%%%%%%%%%%% H ROT %%%%%%%%%%%%%%%%%%%%%%%%%

I2 = rgb2hsv(I); 
H = I2(:,:,1); %% HUE
hue1 = H.*360;

for MEAN_HUE_DIFF = 1
        J1 = rgb2hsv(J1); 
        H = J1(:,:,1); %% HUE
        hue2 = H.*360;

        mn_dif = mean(mean(nonzeros(hue2 - hue1)))
        hue1 = hue1.*double(mk);
        mn_org = mean(mean(nonzeros(hue1)))
end

for OPT_DEG_SEARCH = 1
        deg = 20;
        minMSE = inf;
        optDEG = 0;

        while deg <=360
            break;

            J2 = H_rot(I, mk, deg);

            %msE = immse(uint8(J1), uint8(J2));
            msE = sum(abs(uint8(J1(:)) - uint8(J2(:))));
            if msE < minMSE
                minMSE = msE;
                optDEG = deg; 
            end


            deg = deg + 0.001;

        end
end

optDEG = mn_org + mn_dif %243.501
%minMSE

while 1
    deg = input('Please enter the sigma value for Gaussian Filter: ');
    if deg == -1
        deg = optDEG;
    elseif deg == 0
        break;
    end
    
    figure(3);
    J2 = H_rot(I, mk, deg);
    
    montage([I J2], 'Size', [1 1]);
end





function J = swap_RB(I, mk)
    %%%%%%% SPLIT CHANNELS %%%%%%%%%%%%%%
    R = I(:,:,1);
    G = I(:,:,2);
    B = I(:,:,3);
    
%     figure;
%     montage([R G B], 'Size', [1 1]);
    
    %%%%%%%%%%% MASK CHANNELS %%%%%%%%%%%%
    Rmask = R.*mk;
    Gmask = G.*mk;
    Bmask = B.*mk;

    %%%%%%%%% SWAP R - B %%%%%%%%%%%%%%%
    t = Rmask;
    Rmask = Bmask;
    Bmask = t;

    %%%%%%%%%%%% merge %%%%%%%%%%%%%%%
    mk2 = uint8(~mk);
    R = (R.*mk2) + Rmask;
    G = (G.*mk2) + Gmask;
    B = (B.*mk2) + Bmask;

    J = cat(3, R, G, B);
end


function J = H_rot(I, mk, deg) %%deg degree of rotation
    
    I = rgb2hsv(I);
    
    H = I(:,:,1); %% HUE
    S = I(:,:,2); 
    V = I(:,:,3); 

    %%%%%%%%% Rotate %%%%%%%%%%%%%%%
    hue = H.*360;
    hue = hue + double(deg*mk); %%% rotate only if mk cell is 1
    H = mod(hue, 360) / 360;  %%% taking mod to compensate full rotations
    

    %%%%%%%%%%%% merge %%%%%%%%%%%%%%%
    J = cat(3, H, S, V);
    J = hsv2rgb(J);
%     X = J(:,:,3);
%     min(X(:))
    J = J.*(255);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%