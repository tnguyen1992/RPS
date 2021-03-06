%%%%%%%%% WTC %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
clear all 

srcPath = 'C:\Users\Trinh Nguyen\Documents\MATLAB\RPS\hmrData\Data_FP\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\gcData\Data_FP\';                  % processed data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_sub1_FP.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  numOfPart       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_sub1_FP.mat'));
  end


%% 
for id = numOfPart
  
  % load preprocessed data
  filename_sub1    = sprintf(['RPS_%02d_sub1_FP'], id);
  filename_sub2    = sprintf(['RPS_%02d_sub2_FP'], id);

  fprintf('Load preprocessed data...\n');
    file_path_sub1 = strcat(srcPath, filename_sub1,'.mat');
    file_path_sub2 = strcat(srcPath, filename_sub2,'.mat');

    data_sub1=load(file_path_sub1); 
    data_sub2=load(file_path_sub2);
  %%
% -------------------------------------------------------------------------
% Check order
% -------------------------------------------------------------------------
% Column 7 indicates rest
try
evtRest     = find(data_sub1.s(:, 1) > 0);
end

% Column 9 indicates the beginning of each trial
try
evtTrial  = find(data_sub1.s(:, 3) > 0);
end

% dyad 18 is different as button presses are missing, so rest is 1 and
% trial is 3

if id==18
    evtRest     = find(data_sub1.s(:, 1) > 0);
    evtTrial  = find(data_sub1.s(:, 3) > 0);
end


if size(evtTrial,1)<60
    fprintf('Trial number is different than expected 60!\n');
    missingtrials='yes';
else
    fprintf('Trial number is correct!\n');
    missingtrials='no';
end

% each trial lasts 8 seconds, rest is 60 seconds
    durTrial  = round(8 * ...                % duration task condition
                                  data_sub1.fs - 1);
    durRest     = round(60 * ...                 % duration rest condition
                                  data_sub1.fs - 1);




%% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbo, 2);
hboSub1=data_sub1.hbo;
hboSub2=data_sub2.hbo;

% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hboSub1(1, i)) && ~isnan(hboSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    alpha=0.05;
    max_lag=310;
    x=hboSub1(:,i);
    y=hboSub2(:,i);
    % Does x cause y?
    [F_xy,c_v_xy] = granger_cause(x,y,alpha,max_lag);               
    % Does y cause x?
    [F_yx,c_v_yx] = granger_cause(y,x,alpha,max_lag);               

   
      


granger(i, 1:4) = [F_xy,c_v_xy,F_yx,c_v_yx];

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_FP'], id), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'granger');
  fprintf('Data stored!\n\n');
  clear granger

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PD 

clear all 

srcPath = 'C:\Users\Trinh Nguyen\Documents\MATLAB\RPS\hmrData\Data_PD\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\gcData\Data_PD\';                  % processed data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_sub1_PD.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  numOfPart       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_sub1_PD.mat'));
  end


%% 
for id = numOfPart
  
  % load preprocessed data
  filename_sub1    = sprintf(['RPS_%02d_sub1_PD'], id);
  filename_sub2    = sprintf(['RPS_%02d_sub2_PD'], id);

  fprintf('Load preprocessed data...\n');
    file_path_sub1 = strcat(srcPath, filename_sub1,'.mat');
    file_path_sub2 = strcat(srcPath, filename_sub2,'.mat');

    data_sub1=load(file_path_sub1); 
    data_sub2=load(file_path_sub2);
  %%
% -------------------------------------------------------------------------
% Check order
% -------------------------------------------------------------------------
% Column 7 indicates rest
try
evtRest     = find(data_sub1.s(:, 7) > 0);
end

% Column 9 indicates the beginning of each trial
try
evtTrial  = find(data_sub1.s(:, 9) > 0);
end

% dyad 18 is different as button presses are missing, so rest is 1 and
% trial is 3

if id==18
    evtRest     = find(data_sub1.s(:, 1) > 0);
    evtTrial  = find(data_sub1.s(:, 3) > 0);
end


if size(evtTrial,1)<60
    fprintf('Trial number is different than expected 60!\n');
    missingtrials='yes';
else
    fprintf('Trial number is correct!\n');
    missingtrials='no';
end

% each trial lasts 8 seconds, rest is 60 seconds
    durTrial  = round(8 * ...                % duration task condition
                                  data_sub1.fs - 1);
    durRest     = round(60 * ...                 % duration rest condition
                                  data_sub1.fs - 1);




%% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbo, 2);
hboSub1=data_sub1.hbo;
hboSub2=data_sub2.hbo;
% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hboSub1(1, i)) && ~isnan(hboSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    alpha=0.05;
    max_lag=310;
    x=hboSub1(:,i);
    y=hboSub2(:,i);
    % Does x cause y?
    [F_xy,c_v_xy] = granger_cause(x,y,alpha,max_lag);               
    % Does y cause x?
    [F_yx,c_v_yx] = granger_cause(y,x,alpha,max_lag);               

   
      


granger(i, 1:4) = [F_xy,c_v_xy,F_yx,c_v_yx];

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_PD'], id), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'granger');
  fprintf('Data stored!\n\n');
  clear granger

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PS 

clear all 

srcPath = 'C:\Users\Trinh Nguyen\Documents\MATLAB\RPS\hmrData\Data_PS\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\gcData\Data_PS\';                  % processed data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_sub1_PS.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  numOfPart       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_sub1_PS.mat'));
  end


%% 
for id = numOfPart
  
  % load preprocessed data
  filename_sub1    = sprintf(['RPS_%02d_sub1_PS'], id);
  filename_sub2    = sprintf(['RPS_%02d_sub2_PS'], id);

  fprintf('Load preprocessed data...\n');
    file_path_sub1 = strcat(srcPath, filename_sub1,'.mat');
    file_path_sub2 = strcat(srcPath, filename_sub2,'.mat');

    data_sub1=load(file_path_sub1); 
    data_sub2=load(file_path_sub2);
  %%
% -------------------------------------------------------------------------
% Check order
% -------------------------------------------------------------------------
% Column 7 indicates rest
% try
% evtRest     = find(data_sub1.s(:, 1) > 0);
% end
% 
% % Column 9 indicates the beginning of each trial
% try
% evtTrial  = find(data_sub1.s(:, 3) > 0);
% end
% 
% % dyad 18 is different as button presses are missing, so rest is 1 and
% % trial is 3
% 
% if id==18
%     evtRest     = find(data_sub1.s(:, 1) > 0);
%     evtTrial  = find(data_sub1.s(:, 3) > 0);
% end
% 
% 
% if size(evtTrial,1)<60
%     fprintf('Trial number is different than expected 60!\n');
%     missingtrials='yes';
% else
%     fprintf('Trial number is correct!\n');
%     missingtrials='no';
% end
% 
% % each trial lasts 8 seconds, rest is 60 seconds
%     durTrial  = round(8 * ...                % duration task condition
%                                   data_sub1.fs - 1);
%     durRest     = round(60 * ...                 % duration rest condition
%                                   data_sub1.fs - 1);
% 



%% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbo, 2);
hboSub1=data_sub1.hbo;
hboSub2=data_sub2.hbo;

clear data_sub1 data_sub2

% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hboSub1(1, i)) && ~isnan(hboSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    alpha=0.05;
    max_lag=310;
    x=hboSub1(:,i);
    y=hboSub2(:,i);
    % Does x cause y?
    [F_xy,c_v_xy] = granger_cause(x,y,alpha,max_lag);               
    % Does y cause x?
    [F_yx,c_v_yx] = granger_cause(y,x,alpha,max_lag);               

   
      


granger(i, 1:4) = [F_xy,c_v_xy,F_yx,c_v_yx];

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_PS'], id), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'granger');
  fprintf('Data stored!\n\n');
  clear granger

end