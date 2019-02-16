clc;
clear;
close all;

coords = dlmread('coords.csv');
free_indices = dlmread('free_indices.csv');
no_particles = dlmread('no_particles.csv');
node_info = dlmread('node_info.csv');
properties = dlmread('properties.csv');
X_init = dlmread('X_init.csv');

tic;
[Xk,f_val] = drape_sim(coords,properties,node_info,no_particles,free_indices,X_init);
toc;



% mex -R2018a drape_sim.cpp opt_obj.cpp /usr/local/lib/libnlopt.so