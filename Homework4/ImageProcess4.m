clear;
clc;

test1=imread('Users/macbookpro/Desktop/test1.pgm');
test2=imread('Users/macbookpro/Desktop/test2.tif');

test1Mask3=GaussianFilter(test1,3,1.5);
test1Mask5=GaussianFilter(test1,5,1.5);
test1Mask7=GaussianFilter(test1,7,1.5);

figure
subplot(2,2,1);imshow(test1);title('original');
subplot(2,2,2);imshow(test1Mask3);title('Mask 3*3');
subplot(2,2,3);imshow(test1Mask5);title('Mask 5*5');
subplot(2,2,4);imshow(test1Mask7);title('Mask 7*7');

test2Mask3=GaussianFilter(test2,3,1.5);
test2Mask5=GaussianFilter(test2,5,1.5);
test2Mask7=GaussianFilter(test2,7,1.5);

figure
subplot(2,2,1);imshow(test2);title('original');
subplot(2,2,2);imshow(test2Mask3);title('Mask 3*3');
subplot(2,2,3);imshow(test2Mask5);title('Mask 5*5');
subplot(2,2,4);imshow(test2Mask7);title('Mask 7*7');

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


