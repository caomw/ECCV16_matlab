clear all;
close all;
clc;
mainPath ='/research2/IR_normal_small/save%03d/';
subPath = '%d/';
imPath = '%d.bmp';
normalpath = '12_Normal.bmp';
mean_=0.0;
sum_ = 0.0;
for i = 1:101
    for t =1:9
            normal = imread(sprintf([mainPath,subPath,normalpath],i,t));
            mask = im2bw(normal,0.05) ;
            numpixel = nnz(mask);
        for m=1:12
            im = im2double(imread(sprintf([mainPath,subPath,imPath],i,t,m)));
            im = im*2.0 -1.0;
            im = im.*double(mask);
            sum_ = sum_+ sum(im(:))/numpixel;
        end
    end
end
mean_ = sum_/(101*9*12)