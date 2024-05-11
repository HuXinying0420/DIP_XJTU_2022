clear
clc

test1=imread('C:\Users\LENOVO\Desktop\ImageProcessHomework\Homework5\test1.pgm');
test2=imread('C:\Users\LENOVO\Desktop\ImageProcessHomework\Homework5\test2.tif');
test1 =mat2gray(test1);
test2 =mat2gray(test2);
%Task1
Butterworth_low(test1,200,5)
Butterworth_low(test2,200,5)
Gaussian_low(test1,100)
Gaussian_low(test2,100)

function p=Butterworth_low(Img_in,D0,n)
[M,N]=size(Img_in);    M=2*M;N=2*N;
u = -M/2:(M/2-1);   v = -N/2:(N/2-1);
[u,v] = meshgrid(u,v);
D = sqrt(u.^2+v.^2);
%设计滤波器：
H = 1./(1+(D./D0).^(2*n)); 
Img_fft= fft2(Img_in,size(H,1),size(H,2));
Img_fft_shift = fftshift(Img_fft);
Img_Butterworth = Img_fft_shift.*H;
Img_out = ifft2(ifftshift(Img_Butterworth));
Img_out = Img_out(1:size(Img_in,1),1:size(Img_in,2));
%Img_out =mat2gray(Img_out);
%输出图像：
figure;
subplot(2,2,1);imshow(Img_in);title('原始图像');
subplot(2,2,2);imshow(log(1+abs(Img_fft_shift)),[]);title('频域图像');
subplot(2,2,3);plot3(u,v,H);title('巴特沃斯滤波器');
subplot(2,2,4);imshow(Img_out);title(['滤波后的图像,D0=',num2str(D0),',n=',num2str(n)]);
%计算功率谱比：
S = 0;S1 = 0;
[P,Q] = size(Img_fft_shift);
for a = 1:P
    for b=1:Q
        S1 = S1+(abs(Img_Butterworth(a,b)))^2;
        Img_out = (abs(Img_fft_shift(a,b)))^2;
        S=S+Img_out;
    end
end
p = S1/S;
end

function p=Gaussian_low(Img_in,D0)
[M,N]=size(Img_in);    M=2*M;N=2*N;
u = -M/2:(M/2-1);   v = -N/2:(N/2-1);
[u,v] = meshgrid(u,v);
D = sqrt(u.^2+v.^2);
%设计滤波器：
H = exp(-D.^2/(2.*(D0.^2))); 
Img_fft= fft2(Img_in,size(H,1),size(H,2));
Img_fft_shift = fftshift(Img_fft);
Img_Gaussian = Img_fft_shift.*H;
Img_out = ifft2(ifftshift(Img_Gaussian));
Img_out = Img_out(1:size(Img_in,1),1:size(Img_in,2));
%Img_out =mat2gray(Img_out);
%输出图像：
figure;
subplot(2,2,1);imshow(Img_in);title('原始图像');
subplot(2,2,2);imshow(log(1+abs(Img_fft_shift)),[]);title('频域图像');
subplot(2,2,3);plot3(u,v,H);title('高斯滤波器');
subplot(2,2,4);imshow(Img_out);title(['滤波后的图像,D0=',num2str(D0)]);
%计算功率谱比：
S = 0;S1 = 0;
[P,Q] = size(Img_fft_shift);
for a = 1:P
    for b=1:Q
        S1 = S1+(abs(Img_Gaussian(a,b)))^2;
        Img_out = (abs(Img_fft_shift(a,b)))^2;
        S=S+Img_out;
    end
end
p = S1/S;
end



