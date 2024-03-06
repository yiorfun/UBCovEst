%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 3 AND B2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIGURE 3 Proteomics Dataset %%%
load('Proteomics_FOUR_MATS.mat');
figure;imagesc(PD_107_EST);colormap jet;colorbar; %%% in Figure 3 A third

%%% FIGURE 3 echo-planar spectroscopic imaging Dataset %%%
load('EPSI_FOUR_MATS.mat');
figure;imagesc(EPSI_227_EST);colormap jet;colorbar; %%% in Figure 3 A third

%%% FIGURE B2 %%%

load('PD_EPSI_RESULT.mat')

figure;imagesc(PD_Sigma_samp);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(PD_Sigma_prop);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(PD_Sigma_soft);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(PD_Sigma_full);colormap jet;colorbar;snapnow;caxis([-1,1]);

figure;imagesc(EPSI_Sigma_samp);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(EPSI_Sigma_prop);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(EPSI_Sigma_soft);colormap jet;colorbar;snapnow;caxis([-1,1]);
figure;imagesc(EPSI_Sigma_full);colormap jet;colorbar;snapnow;caxis([-1,1]);

