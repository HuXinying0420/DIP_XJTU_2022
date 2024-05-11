clc
clear

lena=imread('C:\Users\LENOVO\Desktop\ImageProcessHomework\Homework6\lena.bmp');

%task1
lena_g1=GaussianNoise(lena,0,0.1)
lena_g2=GaussianNoise(lena,0,0.3)
lena_g3=GaussianNoise(lena,0,0.5)
lena_g4=GaussianNoise(lena,0.3,0.1)
lena_g5=GaussianNoise(lena,-0.3,0.1)

figure
subplot(2,3,1);imshow(lena);title('原始图像');
subplot(2,3,4);imshow(lena_g1);title('高斯噪声，av=0,std=0.1');
subplot(2,3,5);imshow(lena_g2);title('高斯噪声，av=0,std=0.3');
subplot(2,3,6);imshow(lena_g3);title('高斯噪声，av=0,std=0.5');
subplot(2,3,2);imshow(lena_g4);title('高斯噪声，av=0.3,std=0.1');
subplot(2,3,3);imshow(lena_g5);title('高斯噪声，av=-0.3,std=0.1');

lena_gg=GaussianFilter(lena_g1,5,1.5)
lena_gm=MedianFilter(lena_g1,5)

figure
subplot(2,2,1);imshow(lena);title('原始图像');
subplot(2,2,2);imshow(lena_g1);title('高斯噪声，av=0,std=0.1');
subplot(2,2,3);imshow(lena_gg);title('高斯滤波，5*5，sigma=1.5');
subplot(2,2,4);imshow(lena_gm);title('中值滤波，5*5');

%task2
lena_s1=SaltPepperNoise(lena,0,0.1)
lena_s2=SaltPepperNoise(lena,0.1,0)
lena_s3=SaltPepperNoise(lena,0.1,0.1)

figure
subplot(2,2,1);imshow(lena);title('原始图像');
subplot(2,2,2);imshow(lena_s1);title('椒盐噪声，a=0,b=0.1');
subplot(2,2,3);imshow(lena_s2);title('椒盐噪声，a=0.1,b=0');
subplot(2,2,4);imshow(lena_s3);title('椒盐噪声，a=0.1,b=0.1');

lena_sg=GaussianFilter(lena_s3,5,1.5)
lena_sm=MedianFilter(lena_s3,5)

figure
subplot(2,2,1);imshow(lena);title('原始图像');
subplot(2,2,2);imshow(lena_s3);title('椒盐噪声，a=0.1，b=0.1');
subplot(2,2,3);imshow(lena_sg);title('高斯滤波，5*5，sigma=1.5');
subplot(2,2,4);imshow(lena_sm);title('中值滤波，5*5');

lena_wc1=Contraharmonic(lena_s1,1)
lena_wc2=Contraharmonic(lena_s1,-1)
lena_bc1=Contraharmonic(lena_s2,1)
lena_bc2=Contraharmonic(lena_s2,-1)
lena_sc1=Contraharmonic(lena_s3,1)
lena_sc2=Contraharmonic(lena_s3,-1)
figure
subplot(3,3,1);imshow(lena_s1);title('白盐噪声');
subplot(3,3,2);imshow(lena_wc1);title('白盐噪声，Q=1');
subplot(3,3,3);imshow(lena_wc2);title('白盐噪声，Q=-1');
subplot(3,3,4);imshow(lena_s2);title('胡椒噪声');
subplot(3,3,5);imshow(lena_bc1);title('胡椒噪声，Q=1');
subplot(3,3,6);imshow(lena_bc2);title('胡椒噪声，Q=-1');
subplot(3,3,7);imshow(lena_s3);title('椒盐噪声');
subplot(3,3,8);imshow(lena_sc1);title('椒盐噪声，Q=1');
subplot(3,3,9);imshow(lena_sc2);title('椒盐噪声，Q=-1');

%加高斯噪声：
function Img_out=GaussianNoise(Img,av,std)
[M,N]=size(Img);
u1=rand(M,N);   u2=rand(M,N);
x=std*sqrt(-2*log(u1)).*cos(2*pi*u2)+av;
Img_out=uint8(255*(double(Img)/255+x));
end

%加椒盐噪声：
function Img_out=SaltPepperNoise(Img,a,b)
[M,N]=size(Img);
x=rand(M,N);
Img_out=Img;
Img_out(find(x<=a))=0;
Img_out(find(x>a&x<(a+b)))=255;
end

% 高斯滤波：
function Img_out=GaussianFilter(Img,masksize,sigma)
for i=1:masksize
    for j=1:masksize
        x=i-ceil(masksize/2);
        y=j-ceil(masksize/2);
        h(i,j)=exp(-(x^2+y^2)/(2*sigma^2))/(2*pi*sigma^2);
    end
end
Img_out=uint8(conv2(Img,h,'same'));
end

% 中值滤波：
function Img_out=MedianFilter(Img,masksize)
exsize=floor(masksize/2);   %各方向扩展大小
Imgex=padarray(Img,[exsize,exsize],'replicate','both'); %扩展图片
[m,n]=size(Img);
Img_out=Img;    %将Img_out准备为和Img相同的size
for i=1:m
    for j=1:n
        neighbor=Imgex(i:i+masksize-1,j:j+masksize-1);  %截取邻域
        Img_out(i,j)=median(neighbor(:));   %中值滤波
    end
end
end

%反谐波均值滤波：
function Img_out=Contraharmonic(Img,Q)
[M,N]=size(Img);
ImgSize=3;   ImgSize=(ImgSize-1)/2;
Img_out=Img;
for x=1+ImgSize:1:M-ImgSize
    for y=1+ImgSize:1:M-ImgSize
        is=Img(x-ImgSize:1:x+ImgSize,y-ImgSize:1:y+ImgSize);
        Img_out(x,y)=sum(double(is(:)).^(Q+1))/sum(double(is(:)).^(Q));
    end
end
end
