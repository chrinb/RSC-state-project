function calculate_ca_activity(sData)

% Written by Christoffer Berge | Vervaeke lab

% Compute binazired deconvolved Ca2+ signal based on method in Grosmark et
% al. (2021).

% Load data
dff         = sData.imdata.roiSignals(2).newdff;
deconvolved = sData.imdata.roiSignals(2).ciaDeconvolved;
denoised    = sData.imdata.roiSignals(2).ciaDenoised;

% Compute deconvolution noise: median absolute deviation of the residual  
% of observed trace (dff) and denoised trace reconstruction (denoised)

% loop over ROIs and compute MAD for each ROI
MAD = zeros( size( dff,1), 1);
for roinr = 1:size(dff,1)
    
    observed_trace      = dff(roinr,:);
    denoise_trace_recon = denoised( roinr,:);

    MAD(roinr) = median( abs( observed_trace - denoise_trace_recon ));

end

% Normalize spike estimates (deconvolved trace) by deconvolution noise and
% binarize

mad_threshold = 1.25;

normalized_deconv = zeros(size(deconvolved));
threshold_deconv  = zeros( size(deconvolved));

for roinr = 1:size(deconvolved,1)
    normalized_deconv(roinr,:) = deconvolved(roinr,:) ./ MAD(roinr);

%     threshold_deconv(roinr,:) = 

end



