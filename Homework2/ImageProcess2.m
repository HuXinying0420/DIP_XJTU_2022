clear;
clc;
ImageA=imread('Users/macbookpro/Desktop/Image A.jpg');
ImageB=imread('Users/macbookpro/Desktop/Image B.jpg');
subplot(1,2,1);imshow(ImageA);title('Image A');
subplot(1,2,2);imshow(ImageB);title('Image B');
cpselect(ImageA,ImageB);

T=cp2tform(fixedPoints,movingPoints,'affine');
ImageT=imtransform(ImageB,T);
figure(2);
subplot(1,2,1);imshow(ImageA);title('Image A');
subplot(1,2,2);imshow(ImageT);title('Image T');

