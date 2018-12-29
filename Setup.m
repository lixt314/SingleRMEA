% Run all settings and save figures
clear all
clc
addpath('tools')

%% generate colors and symbols
cols_symbs
%% classify cell types and generate figures
celltype_classification('donor_abc_12k','merged',1);


