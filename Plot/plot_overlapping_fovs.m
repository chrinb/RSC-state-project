function plot_overlapping_fovs(fov1, fov2)

% Written by Christoffer Berge | Vervaeke lab

% plot overlapping FOVs
moving = fov1;
fixed  = fov2;

tformEstimate = imregcorr(moving,fixed, 'similarity');
Rfixed        = imref2d(size(fixed));
movingReg     = imwarp(moving,tformEstimate,"OutputView",Rfixed);

figure, imshowpair(fixed,movingReg,"montage")
figure, imshowpair(fixed,movingReg,"falsecolor");