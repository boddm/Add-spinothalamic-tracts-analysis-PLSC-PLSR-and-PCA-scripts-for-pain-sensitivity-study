clc,clear

%%
load('../Data/mean_dti.mat');


%%
load('data_副本.mat')

brain_data = [mean_dti{1,1}.FA, mean_dti{1,2}.FA ...
    mean_dti{1,1}.AD, mean_dti{1,2}.AD ...
    mean_dti{1,1}.MD, mean_dti{1,2}.MD ...
    mean_dti{1,1}.RD, mean_dti{1,2}.RD];

%%
% Ind = [29 43];
% 8号被试行为学数据有问题
Ind = [26 33];

%%
beh_data(:, 5:6) = [];
beh_data(Ind, :) = [];
brain_data(Ind, :) = [];
diagnosis(Ind, :) = [];

%%
save('data.mat', "diagnosis", "brain_data", "beh_data");