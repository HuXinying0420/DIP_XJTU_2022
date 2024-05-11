clc
clear

lena=imread('C:\Users\LENOVO\Desktop\ImageProcessHomework\Homework6\lena.bmp');
%task3

figure(1);
subplot(241);
imshow(uint8(lena)); title('lena');
subplot(242);
[~, ~, lena_m1] = Motion(lena, 0.1, 0.1, 1);
imshow(lena_m1); title('lena motion');
subplot(246);
lena_m2 = imnoise(lena_m1,'gaussian',0,0.01);
imshow(uint8(lena_m2)); title('lena motion gaussian');

subplot(243);
[~, ~, image_IF] = Weina(lena_m1, 0.03);
imshow(image_IF); title('lena weina');
subplot(244);
[image_IF] = Zuixiao(lena, 0.001, 0);
imshow(uint8(image_IF)); title('lena zuixiao');
subplot(247);
[~, ~, image_IF] = Weina(lena_m2, 0.03);
imshow(image_IF); title('lena weina');
subplot(248);
[image_IF] = Zuixiao(lena, 0.001, 0.001);
imshow(uint8(image_IF)); title('lena zuixiao');

%% hw6 Helper functions
%

% Motion function
function [image_F, image_G, image_IF] = Motion(image, a, b, T)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            x = u - P / 2;
            y = v - Q / 2;
            temp = pi * (x * a + y * b);
            if temp == 0
                temp = 1;
            end
            image_H(u,v) = (T / temp) * sin(temp) * exp(- temp * sqrt(-1));
            image_G(u,v) = image_H(u,v)*image_F(u,v);
        end
    end
    image_IF = uint8(abs(ifft2(ifftshift(image_G))));
end

% Weina function
function [image_F, image_G, image_IF] = Weina(image, K)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    a = 0.1; b = 0.1; T = 1;
    image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            x = u - P / 2;
            y = v - Q / 2;
            temp = pi * (x * a + y * b);
            if temp == 0
                temp = 1;
            end
            image_H(u,v) = (T / temp) * sin(temp) * exp(- temp * sqrt(-1));
        end
    end
    for u = 1:P
        for v = 1:Q
            image_G(u,v)=(1/image_H(u,v)*(abs(image_H(u,v)))^2/((abs(image_H(u,v)))^2+K))*image_F(u,v);
        end
    end
    image_If = ifft2(ifftshift(image_G));
    image_IF = 256.*image_If./max(max(image_If));
    image_IF = uint8(real(image_IF));
end

% Zuixiao function
function [image_IF] = Zuixiao(image, gamma, noise_var)
    [hei,wid,~] = size(image);
    % Simulate a motion blur.
    LEN = 50;
    THETA = 45;
    PSF = fspecial('motion', LEN, THETA);
    blurred = imfilter(image, PSF, 'conv', 'circular');
    % Inverse filter
    Pf = psf2otf(PSF,[hei,wid]);
    % Simulate additive noise.
    noise_mean = 0;
    blurred_noisy = imnoise(blurred, 'gaussian', ...
                            noise_mean, noise_var);
    % Try restoration using  Home Made Constrained Least Squares Filtering.
    p = [0 -1 0;-1 4 -1;0 -1 0];
    P = psf2otf(p,[hei,wid]);
    If = fft2(blurred_noisy);
    numerator = conj(Pf);
    denominator = Pf.^2 + gamma*(P.^2);
    image_IF = ifft2( numerator.*If./ denominator );

%     image_F = fftshift(fft2(double(image)));
%     [P, Q] = size(image_F);
%     a = 0.1; b = 0.1; T = 1;
%     image_H = zeros(P,Q); image_G = zeros(P,Q);
%     for u = 1:P
%         for v = 1:Q
%             x = u - P / 2;
%             y = v - Q / 2;
%             temp = pi * (x * a + y * b);
%             if temp == 0
%                 temp = 1;
%             end
%             image_H(u,v) = (T / temp) * sin(temp) * exp(- temp * sqrt(-1));
%         end
%     end
%     p=[0 -1 0;-1 4 -1;0 -1 0];
%     image_P = psf2otf(p,[P,Q]);
%     image_Hc = conj(image_H);
%     for u = 1:P
%         for v = 1:Q
%             image_G(u,v) = (image_Hc(u,v)*image_F(u,v)/(image_H(u,v))^2+gamma*(image_P(u,v)^2));
%         end
%     end
%     image_IF = ifft2(ifftshift(image_G));
%     image_IF = 256.*image_If./max(max(image_If));
%     image_IF = uint8(real(image_IF));
end

% Helper functions
% image = lena_m1;
% K = 0.0205;
% i = 0;
% Kmin = 0.0205;
% summin = 999999999999;
% for i = 1:100
%     image_F = fftshift(fft2(double(image)));
%     image_H = conj(image_F) / (abs(image_F).^2 + K);
%     image_G = image_F .* image_H;
%     image_IF = uint8(255*abs(ifft2(ifftshift(image_G))));
%     
%     K = K + 0.0001;
%     sumed = sum(sum(double(abs(image - image_IF))));
%     if sumed < summin
%         summin = sumed;
%         Kmin = K;
%     end
% end

%% hw4 Helper functions
%

% Helper function
function image_C=Conv(image, maskSize, Filter)
    [image_M, image_N] = size(image);
    image_extend = wextend('2D','zpd',image,floor(maskSize/2));
    image_C = zeros([image_M, image_N]);
    for i = 1:maskSize
        for j = 1:maskSize
            image_C = image_C + double(image_extend(i:(image_M + i - 1), j:(image_N + j - 1))) * Filter(i, j);
        end
    end
end

% Median Filter by maskSize
function image_M = Median(image, maskSize)
    MFM = ones(maskSize) / maskSize^2;
    image_M = Conv(image, maskSize, MFM);
end

% Gaussian Filter by maskSize
function image_G=Gaussian(image, maskSize)
    GFM = zeros(maskSize);
    sigma = 1.5;
    for i = 1:maskSize
        for j = 1:maskSize
            x = i - floor(maskSize/2) - 1;
            y = j - floor(maskSize/2) - 1;
            GFM(i,j) = exp(-(x^2+y^2)/(2*sigma^2))/(2*pi*sigma^2);
        end
    end
    GFM = GFM / sum(sum(GFM));
    image_G = Conv(image, maskSize, GFM);
end

% Unsharp Filter by maskSize
function image_U = Unsharp(image)
    image_blur = Gaussian(image, 3);
    image_unsharp_mask = double(image) - image_blur;
    image_sharpened = double(image) + image_unsharp_mask;
    image_U = image_sharpened;
end

% Sobel Filter by maskSize
function image_S = Sobel(image)
    SFMx = [1,0,-1;2,0,-2;1,0,-1];
    SFMy = [1,2,1;0,0,0;-1,-2,-1];
    image_Sx = Conv(image, 3, SFMx);
    image_Sy = Conv(image, 3, SFMy);
    image_S = abs(image_Sx) + abs(image_Sy);
end

% Laplace Filter by maskSize
function image_L = Laplace(image)
    LFM = [1,1,1;1,-8,1;1,1,1];
    image_L = Conv(image, 3, LFM);
end

% Canny Filter by maskSize
function image_C = Canny(image)
    image_C = edge(image, 'canny');
end

%% hw5 Helper functions
%

% Helper function
function [image_F, image_G, image_IF, image_r] = BLPF(image, bn, D0)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_Ft = 0; image_Gt = 0;
    image_D = zeros(P,Q); image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            image_D(u,v) = sqrt((u-fix(P/2))^2+(v-fix(Q/2))^2);
            image_H(u,v) = 1/(1+(image_D(u,v)/D0)^(2*bn));
            image_G(u,v) = image_H(u,v)*image_F(u,v);
            image_Fd = (abs(image_F(u,v)))^2;
            image_Ft = image_Ft+image_Fd;
            image_Gd = (abs(image_G(u,v)))^2;
            image_Gt = image_Gt+image_Gd;
        end
    end
    image_IF = uint8(real(ifft2(ifftshift(image_G))));
    image_r = image_Gt/image_Ft;
end

% GLPF function
function [image_F, image_G, image_IF, image_r] = GLPF(image, D0)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_Ft = 0; image_Gt = 0;
    image_D = zeros(P,Q); image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            image_D(u,v) = sqrt((u-fix(P/2))^2+(v-fix(Q/2))^2);
            image_H(u,v) = exp(-image_D(u,v)^2/(2*D0^2));
            image_G(u,v) = image_H(u,v)*image_F(u,v);
            image_Fd = (abs(image_F(u,v)))^2;
            image_Ft = image_Ft+image_Fd;
            image_Gd = (abs(image_G(u,v)))^2;
            image_Gt = image_Gt+image_Gd;
        end
    end
    image_IF = uint8(real(ifft2(ifftshift(image_G))));
    image_r = image_Gt/image_Ft;
end

% BHPF function
function [image_F, image_G, image_IF, image_r] = BHPF(image, bn, D0)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_Ft = 0; image_Gt = 0;
    image_D = zeros(P,Q); image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            image_D(u,v) = sqrt((u-fix(P/2))^2+(v-fix(Q/2))^2);
            image_H(u,v) = 1/(1+(D0/image_D(u,v))^(2*bn));
            image_G(u,v) = image_H(u,v)*image_F(u,v);
            image_Fd = (abs(image_F(u,v)))^2;
            image_Ft = image_Ft+image_Fd;
            image_Gd = (abs(image_G(u,v)))^2;
            image_Gt = image_Gt+image_Gd;
        end
    end
    image_IF = uint8(real(ifft2(ifftshift(image_G))));
    image_r = image_Gt/image_Ft;
end

% GHPF function
function [image_F, image_G, image_IF, image_r] = GHPF(image, D0)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_Ft = 0; image_Gt = 0;
    image_D = zeros(P,Q); image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            image_D(u,v) = sqrt((u-fix(P/2))^2+(v-fix(Q/2))^2);
            image_H(u,v) = 1-exp(-image_D(u,v)^2/(2*D0^2));
            image_G(u,v) = image_H(u,v)*image_F(u,v);
            image_Fd = (abs(image_F(u,v)))^2;
            image_Ft = image_Ft+image_Fd;
            image_Gd = (abs(image_G(u,v)))^2;
            image_Gt = image_Gt+image_Gd;
        end
    end
    image_IF = uint8(real(ifft2(ifftshift(image_G))));
    image_r = image_Gt/image_Ft;
end

% LHPF function
function [image_F, image_G, image_IF, image_r] = LHPF(image, c)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_Ft = 0; image_Gt = 0;
    image_D = zeros(P,Q); image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            image_D(u,v) = sqrt((u-fix(P/2))^2+(v-fix(Q/2))^2);
            image_H(u,v) = 1+c*4*pi^2*image_D(u,v)^2; 
            image_G(u,v) = image_H(u,v)*image_F(u,v);
            image_Fd = (abs(image_F(u,v)))^2;
            image_Ft = image_Ft+image_Fd;
            image_Gd = (abs(image_G(u,v)))^2;
            image_Gt = image_Gt+image_Gd;
        end
    end
    image_IF = uint8(real(ifft2(ifftshift(image_G))));
    image_r = image_Gt/image_Ft;
end

% UHPF function
function [image_F, image_G, image_IF, image_r] = UHPF(image, k1, k2, D0)
    image_F = fftshift(fft2(double(image)));
    [P, Q] = size(image_F);
    image_Ft = 0; image_Gt = 0;
    image_D = zeros(P,Q); image_H = zeros(P,Q); image_G = zeros(P,Q);
    for u = 1:P
        for v = 1:Q
            image_D(u,v) = sqrt((u-fix(P/2))^2+(v-fix(Q/2))^2);
            image_H(u,v) = 1-exp(-image_D(u,v)^2/(2*D0^2));
            image_G(u,v) = (k1+k2*image_H(u,v))*image_F(u,v);
            image_Fd = (abs(image_F(u,v)))^2;
            image_Ft = image_Ft+image_Fd;
            image_Gd = (abs(image_G(u,v)))^2;
            image_Gt = image_Gt+image_Gd;
        end
    end
    image_IF = uint8(real(ifft2(ifftshift(image_G))));
    image_r = image_Gt/image_Ft;
end