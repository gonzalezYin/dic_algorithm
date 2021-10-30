% calculate RMSE for simulation data results
clear;
clc;
close all;


% Get disc
load('/home/dulab/Desktop/dic_result/square20/'+"square20_u.mat");
load('/home/dulab/Desktop/dic_result/square20/'+"square20_v.mat");
displacement = cat(3,u,v);


% Get the centers of the simulation data 
gt_folder = "/home/dulab/Downloads/simulationData/";
folder_name = "square20/";
% file_name = fullfile(gt_folder+folder_name, {'hexagon_0001.mat', 'hexagon_0002.mat', 'hexagon_0003.mat', 'hexagon_0004.mat', 'hexagon_0005.mat',... 
%                                              'hexagon_0006.mat', 'hexagon_0007.mat', 'hexagon_0008.mat', 'hexagon_0009.mat', 'hexagon_0010.mat', ...
%                                              'hexagon_0011.mat', 'hexagon_0012.mat', 'hexagon_0013.mat', 'hexagon_0014.mat', 'hexagon_0015.mat', ...
%                                              'hexagon_0016.mat', 'hexagon_0017.mat', 'hexagon_0018.mat', 'hexagon_0019.mat', 'hexagon_0020.mat'});
file_name = fullfile(gt_folder+folder_name, {'square_0001.mat', 'square_0002.mat'});

for i = 1:2
    gt_imds(i) = load(file_name(i));
end

center1 = reshape(cell2mat(gt_imds(1).centers), [2, 36]);
center2 = reshape(cell2mat(gt_imds(2).centers), [2, 36]);

gt_match = [];
for i = 1:30
    if rem(i, 6) == 0
        % 6 12 18 24 30
    else
        match_tmp = cat(3, center1(:, i), center2(:, i+6));
        gt_match = cat(2, gt_match, match_tmp);
    end
end

disc_ref = [];
for i = 1: size(gt_match, 2)
    % Get the coordinates in the reference, and get the u. v
    tmp = displacement(gt_match(1,i,1), gt_match(2,i,1),:);
    disc_ref = cat(2, disc_ref, tmp);
end

% Get (x+u, y+v)
disc_ref = reshape(disc_ref, [2,25]);
gt_match(:,:,1) = gt_match(:,:,1) + disc_ref(:,:);

% Get RMSE errors
RMSE_errors = [];
% gt_match1 = gt_match(:,:,1).';
% gt_match2 = gt_match(:,:,2).';
% gt_match = cat(3, gt_match1, gt_match2);
% tform_gt = estimateGeometricTransform(gt_match(:,:,1), gt_match(:,:,2), 'affine');
RMSE_errors = mean(sqrt(sum((gt_match(:,:,1)-gt_match(:,:,2)).^2, 2)));







