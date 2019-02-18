clc;
clear;
close all;
rosshutdown;
rosinit;

testclient = rossvcclient('opt_w_cuda');
request = rosmessage(testclient);
request.FileNum = '2';
for i=1:2
tic;
response = call(testclient,request,'Timeout',30);
toc;
a = response.Error
end