clc;
clear;
load('fmriEMCINC.mat')% data
BOLD=fmriMCINC(1:137);
brainNetSet=MVND(BOLD,70,2,lab);%%function MVND