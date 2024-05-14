# DigitalImageProcess_XJTU_2022
# DIP_XJTU_2022

This repository contains assignments from the Digital Image and Video Processing course at Xi'an Jiaotong University for the year 2022. Each assignment covers different techniques and algorithms in image processing.

## Course Description

The Digital Image and Video Processing course provides a comprehensive understanding of various image processing techniques, including image transformation, filtering, enhancement, and noise reduction. The course consists of six assignments that give hands-on experience with these techniques.

## Assignments

### hw1: Basic Image Operations

1. **Introduction to Bmp Image Format**: Explain using 7.bmp as an example.
2. **Gray Level Reduction**: Reduce the gray levels of a 512x512 Lena image from 8 to 1 and display.
3. **Mean and Variance Calculation**: Calculate the mean and variance of the Lena image.
4. **Image Scaling**: Scale the Lena image to 2048x2048 using nearest-neighbor, bilinear, and bicubic interpolation methods.
5. **Image Transformation**: Apply horizontal shear (parameter can be set to 1.5 or chosen by yourself) and 30-degree rotation to Lena and Elain images, then scale to 2048x2048 using nearest-neighbor, bilinear, and bicubic interpolation methods.

### hw2: Image Registration

Given two images (Image A and Image B), randomly select 7 points in each image, calculate the transformation matrix H between the two images, and output the transformed image.

### hw3: Histogram-based Image Enhancement

1. **Histogram Plotting**: Plot the histograms of the provided images.
2. **Histogram Equalization**: Perform histogram equalization on all images, compare the equalized images with the original ones, and analyze the improvements.
3. **Histogram Matching**: Based on the observation of the source image histograms, perform histogram matching to enhance the images.
4. **Local Histogram Enhancement**: Perform 7x7 local histogram enhancement on Elain and Lena images.
5. **Image Segmentation**: Segment Elain and Woman images using histograms.

### hw4: Spatial Domain Filtering

1. **Spatial Low-pass Filtering**: Smooth test images test1 and test2 using Gaussian and median filters with kernel sizes of 3x3, 5x5, and 7x7. Analyze the advantages and disadvantages of each.
2. **Gaussian Filter Generation**: Generate a Gaussian filter with a fixed variance (sigma=1.5) and analyze its advantages and disadvantages.
3. **High-pass Filtering**: Apply unsharp masking, Sobel edge detection, Laplace edge detection, and Canny algorithm to test images test3 and test4. Analyze the advantages and disadvantages of each.

### hw5: Frequency Domain Filtering

1. **Frequency Domain Low-pass Filtering**: Design Butterworth and Gaussian low-pass filters (choose appropriate radii, calculate power spectrum ratio) to smooth test images test1 and test2. Analyze the advantages and disadvantages of each.
2. **Frequency Domain High-pass Filtering**: Design Butterworth and Gaussian high-pass filters to enhance edges in the frequency domain. Choose radii and calculate power spectrum ratio, then apply to test images test3 and test4. Analyze the advantages and disadvantages of each.
3. **Other High-pass Filters**: Apply Laplace and Unmask filters to test images test3 and test4. Analyze the advantages and disadvantages of each.
4. **Comparison of Spatial and Frequency Domain Filtering**: Compare and discuss the relationship between spatial low-pass/high-pass filtering and frequency domain low-pass/high-pass filtering. Analyze whether the results are equivalent in the spatial and frequency domains.

### hw6: Image Noise and Denoising

1. **Gaussian Noise Addition and Denoising**: Add Gaussian noise to the Lena image (specify mean and variance) and restore the image using various filters. Analyze the advantages and disadvantages of each.
2. **Salt-and-Pepper Noise Addition and Denoising**: Add salt-and-pepper noise to the Lena image (density of 0.1) and restore the image using learned filters. Analyze the effect of the contraharmonic filter for Q greater than 0 and less than 0.
3. **Wiener Filter Derivation and Implementation**:
   - Implement the blur filter as per Eq. (5.6-11).
   - Blur the Lena image in the 45-degree direction, T=1.
   - Add Gaussian noise to the blurred Lena image (mean=0, variance=10 pixels) to create a blurred noisy image.
   - Restore the image using Eq. (5.8-6) and (5.9-4) respectively, and analyze the advantages and disadvantages of the algorithms.

<!--


课程作业：图像处理基础算法与各类滤波器的Matlab实现
hw1:
1、Bmp图像格式简介,以7.bmp为例说明；
2、把lena 512*512图像灰度级逐级递减8-1显示；
3、计算lena图像的均值方差；
4、把lena图像用近邻、双线性和双三次插值法zoom到2048*2048；
5、把lena和elain图像分别进行水平shear（参数可设置为1.5，或者自行选择）和旋转30度，并采用用近邻、双线性和双三次插值法zoom到2048*2048；


hw2:图像配准
题目要求：
  要求根据已给的两幅图像，在各幅图像中随机找出7个点，计算出两幅图像之间的转换矩阵H，并且输出转换之后的图像。
注：已给图像分别为Image A和Image B。

hw3:直方图图像增强
共8幅（有数字标号）经变亮或者变暗处理的源图像；
要求：
1.把附件图像的直方图画出； 
2.把所有图像进行直方图均衡；输出均衡后的图像和源图像进行比对；分析改善内容；
3.进一步把图像按照对源图像直方图的观察，各自自行指定不同源图像的直方图，进行直方图匹配，进行图像增强；
4.对elain和lena图像进行7*7的局部直方图增强；
5.利用直方图对图像elain和woman进行分割；
另提供了3幅原始图像：没有数字标号的图像。

hw4
1.空域低通滤波器：分别用高斯滤波器和中值滤波器去平滑测试图像test1和2，模板大小分别是3x3 ， 5x5 ，7x7； 分析各自优缺点；
2.-利用固定方差 sigma=1.5产生高斯滤波器. 附件有产生高斯滤波器的方法； 分析各自优缺点；
3.利用高通滤波器滤波测试图像test3,4：包括unsharp masking, Sobel edge detector, and Laplace edge detection；Canny algorithm.分析各自优缺点；

hw5
1频域低通滤波器：设计低通滤波器包括 butterworth and Gaussian (选择合适的半径，计算功率谱比),平滑测试图像test1和2;分析各自优缺点；
2频域高通滤波器：设计高通滤波器包括butterworth and Gaussian，在频域增强边缘。选择半径和计算功率谱比，测试图像test3,4：分析各自优缺点；
3其他高通滤波器：拉普拉斯和Unmask，对测试图像test3,4滤波；分析各自优缺点；
比较并讨论空域低通高通滤波（Project3）与频域低通和高通的关系；试分析高通、低通滤波器在频域和对应的空域滤波结果是否等效。频域滤波结果如何等效在空频域滤波器。
按标准格式提交报告； 


hw6
1.在测试图像上产生高斯噪声lena图-需能指定均值和方差；并用多种滤波器恢复图像，分析各自优缺点；
2.在测试图像lena图加入椒盐噪声（椒和盐噪声密度均是0.1）；用学过的滤波器恢复图像；在使用反谐波分析Q大于0和小于0的作用；
3.推导维纳滤波器并实现下边要求；
(a) 实现模糊滤波器如方程Eq. (5.6-11).
(b) 模糊lena图像：45度方向，T=1；
(c) 再模糊的lena图像中增加高斯噪声，均值= 0 ，方差=10 pixels 以产生模糊图像；
(d)分别利用方程 Eq. (5.8-6)和(5.9-4)，恢复图像；并分析算法的优缺点.

-->

---

This article was originally published on GitHub by HuXinying0420 in the repository DIP_XJTU_2022 at the following address: https://github.com/HuXinying0420/DIP_XJTU_2022
Please note the following when reprinting this content:

1. Source: This article was originally published on GitHub by HuXinying0420.
2. Original URL: https://github.com/HuXinying0420/DIP_XJTU_2022
