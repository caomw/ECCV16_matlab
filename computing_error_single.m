%computing error depends on epochs

clear all;
close all ;

p = 10; % samples
pp =9; % view point 
folder_set = {'011', '016', '021', '022', '033', '036', '038', '053', '059', '092'};

count =1;   


normalpath = '/research2/IR_normal_small/save';
maskpath = '/research2/IR_normal_small/mask';



    
    aloss_total = 0;
    aloss_total2 = 0;
    err_total = 0;
    err_total2 = 0;
    min_aloss = 100.0;
    max_aloss = -100.0;
    good_pixel1 = 10;
    good_pixel2 = 15;
    good_pixel3 = 20;
    A1_total = 0;
    A2_total = 0;
    A3_total = 0;
    
    folder = fullfile(folder_set{7});
    j =5; % tilt angle        
    im1 = im2double(imread(sprintf('%s%s%s%d%s%s',normalpath,folder,'/',j,'/','12_Normal.bmp')));
    im1 = imresize(im1,0.5);
           
    im2 = im2double(imread('~/Dropbox/ECCV_result/preditced_336000.png'));          
    mask = imread(sprintf('%s%s%s%s%d%s',maskpath,'/',folder,'/',j,'/mask.bmp'));
    mask = imresize(mask,0.5);
    n = sum(sum(mask));
    
    im1 = im1.*2.0-1.0;
    im2 = im2.*2.0-1.0;
    

    
    %%%% Mean Error%%%%%
    im1_ = im1.*repmat(mask,1,1,3);
    im2_ = im2.*repmat(mask,1,1,3);
    err = im1_ - im2_;
    err = abs(err);
    err = err.*double(repmat(mask,1,1,3));
    err_mean = sum(sum(sum(err)));
    err_mean = err_mean./n./3;
    
    err_total = err_total+err_mean;
    
    %%%% Median Error%%%%%
    tmp = err(err~=0);
    err_median = median(tmp);
    err_total2 = err_total2 + err_median;
    
    %%%% Angular Error %%%%
    
    r = size(im1,1);    c = size(im1,2);
    im1_ = im1;
    im2_ = im2;
    
    im1_v = reshape(im1_,[r*c,3]);
    im2_v = reshape(im2_,[r*c,3]);
    
    norm1 = sqrt ( im1_v(:,1).^2 + im1_v(:,2).^2+im1_v(:,3).^2 );
    norm2 = sqrt ( im2_v(:,1).^2 + im2_v(:,2).^2+im2_v(:,3).^2 );
    norm1_3 = repmat(norm1,1,3);
    norm2_3 = repmat(norm2,1,3);
            
    im1_vn = im1_v./norm1_3;
    im2_vn = im2_v./norm2_3;
    
    ang = im1_vn'.*im2_vn';
    ang = sum(ang);
    ang_ = reshape(ang',r,c);
    ang_m = ang_.*double(mask);
    %%%% Convert to radian angles
    %imshow(ang_m);
    
    ang_rad = acosd(ang_m);
    se = strel('line',60,400);
    mask = imerode(mask,se);
    ang_rad_m = ang_rad.*double(mask);
    
    
    ang_rad_m(find(ang_rad_m>30))=30;
    
    aloss = sum(ang_rad_m(:))./sum(mask(:));
    aloss_total = aloss_total+aloss;
    tmp2 = ang_rad_m(ang_rad_m~=0);
    aloss_median = median(tmp2);
    aloss_total2 = aloss_total2 + aloss_median;
    
    A1 = length(find(tmp2<good_pixel1)) ./ length(tmp2) * 100;
    A2 = length(find(tmp2<good_pixel2)) ./ length(tmp2) * 100;
    A3 = length(find(tmp2<good_pixel3)) ./ length(tmp2) * 100;
    
    A1_total = A1_total+A1;
    A2_total = A2_total+A2;
    A3_total = A3_total+A3;
   
          

disp(['Mean of Absolute Error is ' num2str(err_total_mean)]);
disp(['Median of Absolute Error is ' num2str(err_total_median)]);

disp(['Mean of Angular Error is ' num2str(aloss_total_mean)]);
disp(['Median of Angular Error is ' num2str(aloss_total_median)]);
disp(['Max angular error  ' num2str(max_aloss)]);
disp(['Min angular error  ' num2str(min_aloss)]);

disp(['Ratio within 10 deg ' num2str(A1_total_mean)]);
disp(['Ratio within 15 deg ' num2str(A2_total_mean)]);
disp(['Ratio within 20 deg ' num2str(A3_total_mean)]);
disp(['Min angualr loss object ',min_class]);
disp(['Min tilt angular loss object ',min_tilt]);


