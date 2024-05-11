clear
clc
test3=imread('C:\Users\LENOVO\Desktop\ImageProcessHomework\Homework5\test3_corrupt.pgm');
test4=imread('C:\Users\LENOVO\Desktop\ImageProcessHomework\Homework5\test4.tif');

test3 =mat2gray(test3);
test4 =uint8(test4);
test4 =mat2gray(test4);
test4_1=test4(1:512,1:512);
%Task3
Laplacian(test3)
Laplacian(test4_1)
Unmask(test3,20,2)
Unmask(test4_1,20,2)

function Laplacian(Img_in)
[M,N]=size(Img_in);    M=2*M;N=2*N;
u = -M/2:(M/2-1);   v = -N/2:(N/2-1);
[u,v] = meshgrid(u,v);
D = sqrt(u.^2+v.^2);
%����˲�����
H = -4.*pi.^2*(u.^2+v.^2); 
Img_fft= fft2(Img_in,size(H,1),size(H,2));
Img_fft_shift = fftshift(Img_fft);
Img_Laplacian = Img_fft_shift.*H;
Img_out = ifft2(ifftshift(Img_Laplacian));
Img_out = Img_out(1:size(Img_in,1),1:size(Img_in,2));
%���ͼ��
figure;
subplot(2,2,1);imshow(Img_in);title('ԭʼͼ��');
subplot(2,2,2);imshow(log(1+abs(Img_fft_shift)),[]);title('Ƶ��ͼ��');
subplot(2,2,3);plot3(u,v,H);title('������˹�˲���');
subplot(2,2,4);imshow(Img_out);title('�˲����ͼ��');
end

% function Laplacian2(Img_in)
% [M,N]=size(Img_in);    M=2*M;N=2*N;
% u = -M/2:(M/2-1);   v = -N/2:(N/2-1);
% [u,v] = meshgrid(u,v);
% D = sqrt(u.^2+v.^2);
% %����˲�����
% H = -4.*pi.^2*(u.^2+v.^2); 
% Img_fft= fft2(Img_in,size(H,1),size(H,2));
% Img_fft_shift = fftshift(Img_fft);
% img_temp=log(1+abs(Img_fft_shift));
% Img_Laplacian = Img_fft_shift.*H;
% Img_out = ifft2(ifftshift(Img_Laplacian));
% Img_out = Img_out(1:size(Img_in,1),1:size(Img_in,2));
% %���ͼ��
% figure;
% subplot(2,2,1);imshow(Img_in(:,1:512));title('ԭʼͼ��');
% subplot(2,2,2);imshow(img_temp(1:512,1:512),[]);title('Ƶ��ͼ��');
% subplot(2,2,3);plot3(u,v,H);title('������˹�˲���');
% subplot(2,2,4);imshow(Img_out(1:512,1:512));title('�˲����ͼ��');
% end

function Unmask(Img_in,D0,n)
Img_fft_shift=fftshift(fft2(Img_in));
[M,N]=size(Img_fft_shift);
for i=1:M
    for j=1:N
        d=sqrt((i-fix(M/2))^2+(j-fix(N/2))^2);
        if d==0
            H(i,j)=0;
        else
            H(i,j)=1/(1+0.414*(D0/d)^(2*n));
        end
        Img_out(i,j)=(1+H(i,j))*Img_fft_shift(i,j);
    end
end
Img_out=real(ifft2(ifftshift(Img_out)));
%���ͼ��
figure;
subplot(2,2,1);imshow(Img_in);title('ԭʼͼ��');
subplot(2,2,2);imshow(log(1+abs(Img_fft_shift)),[]);title('Ƶ��ͼ��');
subplot(2,2,3);
u=-M/2:(M/2-1);v=-N/2:(N/2-1);
[u,v]=meshgrid(u,v);plot3(u,v,H);title('Unmask�˲���');
subplot(2,2,4);imshow(Img_out);title(['�˲����ͼ��,D0=',num2str(D0),',n=',num2str(n)]);
end

% function Unmask2(Img_in,D0,n)
% Img_fft_shift=fftshift(fft2(Img_in));
% [M,N]=size(Img_fft_shift);
% for i=1:M
%     for j=1:N
%         d=sqrt((i-fix(M/2))^2+(j-fix(N/2))^2);
%         if d==0
%             H(i,j)=0;
%         else
%             H(i,j)=1/(1+0.414*(D0/d)^(2*n));
%         end
%         Img_out(i,j)=(1+H(i,j))*Img_fft_shift(i,j);
%     end
% end
% Img_out=real(ifft2(ifftshift(Img_out)));
% temp=log(1+abs(Img_fft_shift));
% %���ͼ��
% figure;
% subplot(2,2,1);imshow(Img_in(:,1:512));title('ԭʼͼ��');
% subplot(2,2,2);imshow(temp(1:512,1:512),[]);title('Ƶ��ͼ��');
% subplot(2,2,3);
% u=-M/2:(M/2-1);v=-N/2:(N/2-1);
% [u,v]=meshgrid(u,v);plot3(u,v,H);title('Unmask�˲���');
% subplot(2,2,4);imshow(Img_out(1:512,1:512));title(['�˲����ͼ��,D0=',num2str(D0),',n=',num2str(n)]);
% end