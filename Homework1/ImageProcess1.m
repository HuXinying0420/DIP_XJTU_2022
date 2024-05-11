clear;
clc;

%Task1
lena=imread('Users/macbookpro/Desktop/lena.bmp');
lena7=lena/2;   lena6=lena/4;   lena5=lena/8;
lena4=lena/16;  lena3=lena/32;  lena2=lena/64;
lena1=lena/128;
subplot(2,4,1);    imshow(lena,[0,256]);title('原图');
subplot(2,4,2);    imshow(lena7,[0,127]);title('7灰度级');
subplot(2,4,3);    imshow(lena6,[0,63]);title('6灰度级');
subplot(2,4,4);    imshow(lena5,[0,31]);title('5灰度级');
subplot(2,4,5);    imshow(lena4,[0,15]);title('4灰度级');
subplot(2,4,6);    imshow(lena3,[0,7]);title('3灰度级');
subplot(2,4,7);    imshow(lena2,[0,3]);title('2灰度级');
subplot(2,4,8);    imshow(lena1,[0,1]);title('1灰度级');

%Task2
mean = mean2(lena)
std = std2(lena);
var = std^2

%Task3
nearest=imresize(lena,[2048,2048],'nearest');
bilinear=imresize(lena,[2048,2048],'bilinear');
bicubic=imresize(lena,[2048,2048],'bicubic');
imtool(nearest);
imtool(bilinear);
imtool(bicubic);

%Task4
elain=imread('Users/macbookpro/Desktop/elain1.bmp');
%构建仿射矩阵：
T_shear=[1,1.5,0;0,1,0;0,0,1];
T_rotation=[cos(pi/6),sin(pi/6),0;-sin(pi/6),cos(pi/6),0;0,0,1];

shear=affine2d(T_shear);
rotation=affine2d(T_rotation);

lena_shear=imwarp(lena,shear);
elain_shear=imwarp(elain,shear);
lena_rotation=imwarp(lena,rotation);
elain_rotation=imwarp(elain,rotation);

figure(1);imshow(lena_shear);title("lena shear");
figure(2);imshow(elain_shear);title("elain shear");
figure(3);imshow(lena_rotation);title("lena rotation");
figure(4);imshow(elain_rotation);title("elain ratation");

lena_shear_zoom1=imresize(lena_shear,[2048,2048],'nearest');
lena_shear_zoom2=imresize(lena_shear,[2048,2048],'bilinear');
lena_shear_zoom3=imresize(lena_shear,[2048,2048],'bicubic');
elain_shear_zoom1=imresize(elain_shear,[2048,2048],'nearest');
elain_shear_zoom2=imresize(elain_shear,[2048,2048],'bilinear');
elain_shear_zoom3=imresize(elain_shear,[2048,2048],'bicubic');

lena_rotation_zoom1=imresize(lena_rotation,[2048,2048],'nearest');
lena_rotation_zoom2=imresize(lena_rotation,[2048,2048],'bilinear');
lena_rotation_zoom3=imresize(lena_rotation,[2048,2048],'bicubic');
elain_rotation_zoom1=imresize(elain_rotation,[2048,2048],'nearest');
elain_rotation_zoom2=imresize(elain_rotation,[2048,2048],'bilinear');
elain_rotation_zoom3=imresize(elain_rotation,[2048,2048],'bicubic');

imtool(lena_shear_zoom1);
imtool(lena_shear_zoom2);
imtool(lena_shear_zoom3);
imtool(elain_shear_zoom1);
imtool(elain_shear_zoom2);
imtool(elain_shear_zoom3);
imtool(lena_rotation_zoom1);
imtool(lena_rotation_zoom2);
imtool(lena_rotation_zoom3);
imtool(elain_rotation_zoom1);
imtool(elain_rotation_zoom2);
imtool(elain_rotation_zoom3);


