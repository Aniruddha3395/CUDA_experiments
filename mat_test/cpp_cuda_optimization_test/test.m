%****************************************************************************************
%
% Author : Aniruddha Shembekar, University of Southern California
%
%****************************************************************************************

clc;
clear;
close all;


mexcuda -c opt_obj.cu -lcudart; 
mex -L/usr/local/cuda/lib64 -lcudart -I./ -R2018a drape_sim.cpp opt_obj.o /usr/local/lib/libnlopt.so;  

file_num = 3;

coords = dlmread(strcat('coords',num2str(file_num),'.csv'));
free_indices = dlmread(strcat('free_indices',num2str(file_num),'.csv'));
no_particles = dlmread(strcat('no_particles',num2str(file_num),'.csv'));
node_info = dlmread(strcat('node_info',num2str(file_num),'.csv'));
properties = dlmread(strcat('properties',num2str(file_num),'.csv'));
X_init = dlmread(strcat('X_init',num2str(file_num),'.csv'));


tic;
[Xk,f_val] = drape_sim(coords,properties,node_info,no_particles,free_indices,X_init);
toc;

f_val
% MEX with only nlopt
% mex -R2018a drape_sim.cpp opt_obj.cpp /usr/local/lib/libnlopt.so

% MEX with nlopt and cuda together
% mexcuda -c err_fun_cuda.cu -lcudart; 
% mex -L/usr/local/cuda/lib64 -lcudart -I./ -R2018a drape_sim.cpp opt_obj.cpp err_fun_cuda.o /usr/local/lib/libnlopt.so;  

%% NEW
% mexcuda -c opt_obj.cu -lcudart; 
% mex -L/usr/local/cuda/lib64 -lcudart -I./ -R2018a drape_sim.cpp opt_obj.o /usr/local/lib/libnlopt.so;  

