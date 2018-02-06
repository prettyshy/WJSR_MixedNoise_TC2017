clear
clc
Original_image_dir  =    'C:\Users\csjunxu\Desktop\TWSCGIN\cleanimages\';
Sdir = regexp(Original_image_dir, '\', 'split');
fpath = fullfile(Original_image_dir, '*.png');
im_dir  = dir(fpath);
im_num = length(im_dir);
method = 'WJSR';
write_MAT_dir = ['C:/Users/csjunxu/Desktop/TWSCGIN/'];
write_sRGB_dir = [write_MAT_dir method];
if ~isdir(write_sRGB_dir)
    mkdir(write_sRGB_dir)
end
Type = 0;
for nSig  =   [10 20 30]         % The standard variance of the additive Gaussian noise;
    for sp =  [.1 .3 .5]           % salt and pepper
        % record all the results in each iteration
        PSNR = zeros(im_num, 1, 'double');
        SSIM = zeros(im_num, 1, 'double');
        T512 = [];
        T256 = [];
        for i = 1:im_num
            ff2=double( imread(fullfile(Original_image_dir, im_dir(i).name)) );
            S = regexp(im_dir(i).name, '\.', 'split');
            
            IMin0=ff2;
            IMin0=double(IMin0);
            % nSig=10;
            % Nimg0=double(IMin0)+nSig*randn(size(IMin0));
            % sp=0.2;
            % [Nimg1]          =   impulsenoise(Nimg0,sp);
            if Type == 0
                imname = sprintf([write_MAT_dir 'noisyimages/G' num2str(nSig) '_SPIN' num2str(sp) '_' im_dir(i).name]);
                %                                     imwrite(Par.nim/255,imname);
                Nimg1 = double( imread(imname));
            elseif Type == 1
                imname = sprintf([write_MAT_dir 'noisyimages/G' num2str(nSig) '_RVIN' num2str(sp) '_' im_dir(i).name]);
                %                                     imwrite(Par.nim/255,imname);
                Nimg1 = double( imread(imname));
            else
                break;
            end
            
            
            K=256;
            n=8;
            bb=8;
            
            IMin0=double(IMin0);
            Nimg=double(Nimg1);
            % [PSNRnois, MSEnois] = psnr(double(IMin0),Nimg)
            fprintf('%s :\n',im_dir(i).name);
            fprintf('The initial value of PSNR = %2.4f, SSIM = %2.4f \n', csnr( Nimg, IMin0, 0, 0 ), cal_ssim( Nimg, IMin0, 0, 0 ));
            
            [AdpMedDenoised,flagDctNImg]           =   adpmedft(Nimg,19);  % for salt & pepper noise
            flagDctImg=flagDctNImg;
            
            outpDct=Nimg;
            if (nSig<=10)
                MaxIterDCT=7;
            elseif (nSig<=20)
                MaxIterDCT=8;
            elseif (nSig<=30)
                MaxIterDCT=10;
            else
                MaxIterDCT=13;
            end
            for k=1:MaxIterDCT  % 7 for nSig=10, 8 for nSig=20, 10 for nSig=30
                preImgDCT=outpDct;
                [outpDct,TDict]=ImpNoisRemove_DCT(IMin0,outpDct,1-flagDctImg);
                outpDct=max(min(outpDct,255),0);
                curImgDCT=outpDct;
                res1=sqrt(sum(sum((curImgDCT-preImgDCT).^2)));
                res2=sqrt(sum(sum(preImgDCT.^2)));
                dif=res1/res2;
                if (dif<=0.013)  % 0.016 for gaussian=10
                    break;
                end
                % figure(3), imshow(uint8(outpDct));
                % title('denoised image by weighted SR with DCT dictionary')
                [AdpMedDenoised,flagDctImg]           =   adpmedft(curImgDCT,19);  % for salt & pepper noise
            end
            
            %---------the parameter structure----------
            par.step         =   2;
            par.nblk         =   5;
            par.patchsize    =   bb;
            par.const        =   sqrt(1.15);
            par.sig          =   nSig;
            %------------------------------------------
            refImg=outpDct;
            dict=TDict;
            flagNim=1-flagDctNImg;
            preDim=Nimg;
            par.sig=max(5,par.sig/1.4);
            tau=nSig;
            Lambda=sqrt(nSig);
            Lambda2=4;
            if (nSig<=10)
                MaxIter=3;
            elseif (nSig<=20)
                MaxIter=6;
            elseif (nSig<=30)
                MaxIter=9;
            else
                MaxIter=12;
            end
            for k=1:MaxIter
                NLWeightSRiter=k
                
                [outpDimg,outpTemp]=NLGroupWeightSR(preDim,refImg,par,dict,flagNim,Lambda,Lambda2);
                curDim=outpDimg;
                outpDimg=max(min(outpDimg,255),0);
                % [PSNRdim, MSEdim] = psnr(double(IMin0),outpDimg)
                % [missim_dim, ssim_dim] = ssim_index(double(IMin0),outpDimg);
                % missim_dim
                % figure(8), imshow(uint8(outpDimg));
                % title('denoised image by weighted simutaneous SR with DCT dictionary and global regularization')
                
                resd1=sqrt(sum(sum((curDim-preDim).^2)));
                resd2=sqrt(sum(sum(preDim.^2)));
                dif=resd1/resd2;
                if (dif<0.006)
                    break;
                end
                preDim=outpDimg;
                refImg=outpTemp;
                [AdpMedDenoised,flagDctImg]           =   adpmedft(outpDimg,19);
                flagNim=1-flagDctImg;
                par.sig=max(5,par.sig/2.5);
                Lambda=max(sqrt(nSig)/4,Lambda/2);
            end
            %% calculate the PSNR
            PSNR(i, 1)  =   csnr( outpDimg, IMin0, 0, 0 );
            SSIM(i, 1)      =  cal_ssim( outpDimg, IMin0, 0, 0 );
            imname = sprintf([write_sRGB_dir '/' method '_AMF1_GSPIN_nl' num2str(nSig) '_sp' num2str(sp) im_dir(i).name]);
            imwrite(outpDimg/255,imname);
            fprintf('%s : PSNR = %2.4f, SSIM = %2.4f \n', im_dir(i).name, PSNR(i, 1), SSIM(i, 1) );
        end
        mPSNR=mean(PSNR);
        mSSIM=mean(SSIM);
        name = sprintf([write_MAT_dir '/' method '_AMF1_GSPIN_nl' num2str(nSig) '_sp' num2str(sp) '.mat']);
        save(name,'nSig','sp','PSNR','SSIM','mPSNR','mSSIM');
    end
end



