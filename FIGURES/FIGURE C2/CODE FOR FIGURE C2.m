%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE C2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('Supp_C4_Plots.mat')
figure;imagesc(SIGMAU0);colormap jet;colorbar; 
figure;imagesc(SU);colormap jet;colorbar; 
figure;imagesc(SIGMAU_EST);colormap jet;colorbar; 

figure;imagesc(SIGMAH0);colormap jet;colorbar; 
figure;imagesc(SH);colormap jet;colorbar; 
figure;imagesc(SIGMAH_EST);colormap jet;colorbar; 

figure;imagesc(SIGMAI0);colormap jet;colorbar; 
figure;imagesc(SI);colormap jet;colorbar; 
figure;imagesc(SIGMAI_EST);colormap jet;colorbar; 
