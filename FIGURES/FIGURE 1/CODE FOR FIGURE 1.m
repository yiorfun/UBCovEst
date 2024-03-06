%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MATLAB CODE FOR FIGURE 1 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% FIGURE 1 A Spellman's Dataset %%%
load('Spellman_FOUR_MATS.mat');
figure;imagesc(Spell_RAW);colormap jet;colorbar; %%% in Figure 1 A Left
figure;imagesc(Spell_ALL);colormap jet;colorbar; %%% in Figure 1 A Right

%%% FIGURE 1 B Proteomics Dataset %%%
load('Proteomics_FOUR_MATS.mat');
figure;imagesc(PD_RAW);colormap jet;colorbar; %%% in Figure 1 B Left
figure;imagesc(PD_ALL);colormap jet;colorbar; %%% in Figure 1 B Right

%%% FIGURE 1 C Seed Quality Dataset %%%
load('Seed_FOUR_MATS.mat');
figure;imagesc(Seed_RAW);colormap jet;colorbar; %%% in Figure 1 C Left
figure;imagesc(Seed_ALL);colormap jet;colorbar; %%% in Figure 1 C Right

%%% FIGURE 1 D Nuclear Magnetic Resonance Dataset %%%
load('NMR_E_FOUR_MATS.mat');
figure;imagesc(NMR_E_RAW);colormap jet;colorbar; %%% in Figure 1 D Left
figure;imagesc(NMR_E_ALL);colormap jet;colorbar; %%% in Figure 1 D Right

%%% FIGURE 1 E Exposome Dataset %%%
load('Exposome_FOUR_MATS.mat');
figure;imagesc(Exposome_RAW);colormap jet;colorbar; %%% in Figure 1 E Left
figure;imagesc(Exposome_ALL);colormap jet;colorbar; %%% in Figure 1 E Right

%%% FIGURE 1 F Metabolites Dataset %%%
load('Metabolites_FOUR_MATS.mat');
figure;imagesc(Metabolites_RAW);colormap jet;colorbar; %%% in Figure 1 F Left
figure;imagesc(Metabolites_ALL);colormap jet;colorbar; %%% in Figure 1 F Right

%%% FIGURE 1 G echo-planar spectroscopic imaging Dataset %%%
load('EPSI_FOUR_MATS.mat');
figure;imagesc(EPSI_RAW);colormap jet;colorbar; %%% in Figure 1 G Left
figure;imagesc(EPSI_ALL);colormap jet;colorbar; %%% in Figure 1 G Right





