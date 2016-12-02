%computing error withlights

clear all;
close all ;

p = 10; % samples
pp =9; % view point(titt angles)
l=12; % number of lights directions 
folder_set = {'011', '016', '021', '022', '033', '036', '038', '053', '059', '092'};

count =1;
values=[];
objects=[];
view=[];


model ='DCGAN.model-51612';
outputpath ='/research2/ECCV_journal/deconv/NIR_single_dis';
suboutputpath ='Deconv_L1_result';
normalpath = '/research2/IR_normal_small/save';
maskpath = '/research2/IR_normal_small/mask';
fileID = fopen(sprintf('~/Dropbox/ECCV_result/deconv/nodetail_2D_pair.csv'),'w'); %show Total mean etc 
fileID2 = fopen(sprintf('~/Dropbox/ECCV_result/deconv/nodetail_2D_pair.csv'),'w'); %show depends on lights 

fprintf(fileID,'mean_error,median_error,mean_angle_error,10deg\n');
fprintf(fileID2,'mean_error,median_error,mean_angle_error,10deg\n');
lights = dir(fullfile(outputpath,suboutputpath,'011/1',model,'*.bmp'));
lights = sort({lights.name});
lights = natsort(lights);

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
    
abs_mean= zeros(l,1);
abs_median=zeros(l,1);
ang_loss=zeros(l,1);
deg_10=zeros(l,1);

loaded = false;
for e =1:length(lights) % lights
    for i = 1: p %object
        folder = fullfile(folder_set{i});
        for j=1:pp % tilt
            sprintf('Processsing object:%d/%d tilt: %d/%d lights: %d/%d \n',i,p,j,pp,e,length(lights))         
            im1 = im2double(imread(sprintf('%s%s%s%d%s%s',normalpath,folder,'/',j,'/','12_Normal.bmp')));
            im1 = imresize(im1,0.5);
            im1 = im1.*2-1;
            mask = imread(sprintf('%s%s%s%s%d%s',maskpath,'/',folder,'/',j,'/mask.bmp'));
            mask = imresize(mask,0.5);
            n = sum(sum(mask));
          
            im2 = im2double(imread(fullfile(outputpath,sprintf('%s/%s',suboutputpath,folder),int2str(j),model,lights{e})));
            im2 = im2.*2-1;
            %%%%% Mean Error%%%%%
            im1_ = im1.*repmat(mask,1,1,3);
            im2_ = im2.*repmat(mask,1,1,3);
            err = im1_ - im2_;
            err = abs(err);
            err = err.*double(repmat(mask,1,1,3));
            err_mean = sum(sum(sum(err)));
            err_mean = err_mean./n./3;
            err_total = err_total+err_mean;
            abs_mean(e)  =  abs_mean(e)+err_mean;
            %%%%%%%% Median Error%%%%%%%
            tmp = err(err~=0);
            err_median = median(tmp);
            if isnan(err_median)
                disp('NaN');
            end
            err_total2 = err_total2 + err_median;
            abs_median(e)  =  abs_median(e)+err_median;
            %%%%%%% Angular Error %%%%%%%%
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
            ang_rad = acosd(ang_m);
            se = strel('line',60,400);
            mask = imerode(mask,se);
            ang_rad_m = ang_rad.*double(mask);
            ang_rad_m(find(ang_rad_m>30))=30;
            
            
            aloss = sum(ang_rad_m(:))./sum(mask(:));
            aloss_total = aloss_total+aloss;
            ang_loss(e) = ang_loss(e) + aloss;
            
            tmp2 = ang_rad_m(ang_rad_m~=0);
            aloss_median = median(tmp2);
            aloss_total2 = aloss_total2 + aloss_median;
            
            A1 = length(find(tmp2<good_pixel1)) ./ length(tmp2) * 100;
            A2 = length(find(tmp2<good_pixel2)) ./ length(tmp2) * 100;
            A3 = length(find(tmp2<good_pixel3)) ./ length(tmp2) * 100;
            
            A1_total = A1_total+A1;
            A2_total = A2_total+A2;
            A3_total = A3_total+A3;
            
            deg_10(e)  = deg_10(e) + A1;
            
%             if aloss > max_aloss
%                 max_aloss = aloss;
%                 max_class = folder_set{i};
%                 max_tilt = j
%             end
%             
%             if min_aloss > aloss
%                 min_aloss = aloss;
%                 min_class =folder_set{i};
%                 min_tilt = j;
%             end
        end
    end
end

 
 err_total_mean = err_total/(p*pp*l);
 err_total_median = err_total2/(p*pp*l);
 aloss_total_mean = aloss_total/(p*pp*l);
 aloss_total_median = aloss_total2/(p*pp*l);
 A1_total_mean = A1_total/(p*pp*l);
 A2_total_mean = A2_total/(p*pp*l);
 A3_total_mean = A3_total/(p*pp*l);
 
 abs_mean= abs_mean/(p*pp);
 abs_median=abs_median/(p*pp); 
 ang_loss= ang_loss/(p*pp);
 deg_10=deg_10/(p*pp);
 
 fprintf(fileID,'%.6f, %.6f,%.6f\n',err_total_mean,aloss_total_mean,A1_total_mean);
 fclose(fileID);
 
 
 for nn= 1:12 %lights
    fprintf(fileID2,'%.6f, %.6f,%.6f,%.6f\n',abs_mean(nn),abs_median(nn),ang_loss(nn),deg_10(nn));
 end
  fclose(fileID2);
fclose('all');
 
 
% disp(['Mean of Absolute Error is ' num2str(err_total_mean)]);
% disp(['Median of Absolute Error is ' num2str(err_total_median)]);
% 
% disp(['Mean of Angular Error is ' num2str(aloss_total_mean)]);
% disp(['Median of Angular Error is ' num2str(aloss_total_median)]);
% disp(['Max angular error  ' num2str(max_aloss)]);
% disp(['Min angular error  ' num2str(min_aloss)]);
% 
% disp(['Ratio within 10 deg ' num2str(A1_total_mean)]);
% disp(['Ratio within 15 deg ' num2str(A2_total_mean)]);
% disp(['Ratio within 20 deg ' num2str(A3_total_mean)]);



