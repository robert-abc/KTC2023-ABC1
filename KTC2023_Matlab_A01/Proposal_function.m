%This script runs the main.m function for each file in the input folder
close all; 
clear all; 
clc;
path(path,'MiscCodes/')
main('TrainingData', 'Output', 7);




%% Reconstructions assessment, for testing purposes only

% score = 0;
% initial_file = 1;
% final_file = 4;
% for ii = initial_file:final_file
%     load(['Output/' num2str(ii) '.mat']);
%     load(['GroundTruths/true' num2str(ii) '.mat']);
%     s = scoringFunction(truth, reconstruction);
%     disp(['Score from target ' num2str(ii) ' = ' num2str(s)])
%     score = score + s;
% 
%     %Ground truth and reconstruction visualizations
%     figure;
%     subplot(1,2,1); imagesc(truth); axis square; clim([0 2]); colorbar; title('Ground truth')
%     subplot(1,2,2); imagesc(reconstruction); axis square; clim([0 2]); colorbar;title('Reconstruction')
% end
% 
% disp(['Final score: ' num2str(score/(final_file-initial_file+1)) ])

