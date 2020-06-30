load('case4_50x50Result.mat');
load('case4Model.mat');

figure;
pdeplot(model, 'XYData', results.NodalSolution(:,1,1));
colormap('jet');

figure;
pdeplot(model, 'XYData', results.NodalSolution(:,1,5));
colormap('jet');

figure;
pdeplot(model, 'XYData', results.NodalSolution(:,1,10));
colormap('jet');


figure;
pdeplot(model, 'XYData', results.NodalSolution(:,1,end));
colormap('jet');

% figure;
% pdeplot(model, 'XYData', results.NodalSolution(:,1,2));
% colormap('jet');
% 
% figure;
% pdeplot(model, 'XYData', results.NodalSolution(:,1,3));
% colormap('jet');

% [numX, numY, numZ] = size(results.NodalSolution);
% for idx = 1:numZ
%     figure;
%     pdeplot(model, 'XYData', results.NodalSolution(:,1,idx));
%     colormap('jet');
% end
