% % % %%%%%%%%% WTC %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
clear all 

srcPath = 'P:\projects\RPS\RPS\procData\hmrData\Data_C\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_C\';                  % processed data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_sub1_C.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  numOfPart       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_sub1_C.mat'));
  end


%%
numOfPermutations=1000; 
for id = numOfPart
  
  % load preprocessed data
  for n=1:1:numOfPermutations
        fprintf('<strong>Permutation %d</strong>\n', n);
  f=0;
  while(f==0)
  try
  filename_sub1    = sprintf(['RPS_%02d_sub1_C'], id);
  k=randi(numOfSources);
  filename_sub2    = sprintf(['RPS_%02d_sub2_C'], numOfPart(k));

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



%% Calc the period 
hbrSub1=data_sub1.hbr;
hbrSub2=data_sub2.hbr;

poi=[6 20];
pnoi = zeros(2,1);

  sigPart1 = [data_sub1.t, hbrSub1(:,1)];
  sigPart2 = [data_sub2.t, hbrSub2(:,1)];
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > poi(1), 1, 'first');
  pnoi(2) = find(period < poi(2), 1, 'last');
                                                                     

% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbr, 2);
coherences  = zeros(numOfChan, 2);
meanCohControl = zeros(1, length(evtTrial));                         % mean coherence in a defined spectrum for condition collaboration      
meanCohRest   = zeros(1, length(evtRest));                              % mean coherence in a defined spectrum for condition baseline
% mean coherence in a defined spectrum for condition baseline
Rsq{numOfChan} = [];
Rsq(:) = {NaN(length(period), length(data_sub1.t))};

% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hbrSub1(1, i)) && ~isnan(hbrSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    sigPart1 = [data_sub1.t, hbrSub1(:,i)];
    sigPart2 = [data_sub2.t, hbrSub2(:,i)];
    [Rsq{i}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
   
      for j=1:1:length(coi)
        Rsq{i}(period >= coi(j), j) = NaN;
      end
   
    
    % calculate mean activation in frequency band of interest
    % cooperation condition
    for j=1:1:length(evtTrial)
        meanCohControl(j)  = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtTrial(j):evtTrial(j) + ...
                            durTrial)));
    end
 

       % rest condition
    for j=1:1:length(evtRest)
        try
        meanCohRest(j)    = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtRest(j):evtRest(j) + ...
                            durRest)));
        catch
           meanCohRest(j) = NaN; 
        end
    end
    
    control    = nanmean(meanCohControl);                                % average mean coherences over trials
    control_1  = nanmean(meanCohControl(1:size(meanCohControl,2)/2));
    control_2  = nanmean(meanCohControl(size(meanCohControl,2)/2+1:end));
    rest       = nanmean(meanCohRest);


coherences(i, 1:2) = [control,rest];
coherences1(i,1:3) = [control_1,control_2,rest];   

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_C_hbr_p%04d'], id, n), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'coherences','coherences1','missingtrials');
  fprintf('Data stored!\n\n');
  clear coherences
  f=1;
  catch
       
   fprintf("Repeat loop iteration");
   f=0;
  end %end for try/catch 
  end %end for while function
  end %end for permutation loop

end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FP 

clear all 

srcPath = 'P:\projects\RPS\RPS\procData\hmrData\Data_FP\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_FP\';                  % processed data location
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
numOfPermutations=1000; 
for id = numOfPart
  
  % load preprocessed data
   for n=1:1:numOfPermutations
        fprintf('<strong>Permutation %d</strong>\n', n);
  f=0;
  while(f==0)
  try
  filename_sub1    = sprintf(['RPS_%02d_sub1_FP'], id);
  k=randi(numOfSources);
  filename_sub2    = sprintf(['RPS_%02d_sub2_FP'], numOfPart(k));

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
evtRest     = find(data_sub1.s(:, 1) > 0);

% Column 9 indicates the beginning of each trial
evtTrial  = find(data_sub1.s(:, 3) > 0);
if size(evtTrial,1)<60
    fprintf('Trial number is different than expected 60! \n');
    missingtrials='yes';
else
    fprintf('Trial number is correct! \n');
    missingtrials='no';
end

% each trial lasts 9 seconds, rest is 60 seconds
    durTrial  = round(8 * ...                % duration task condition
                                  data_sub1.fs - 1);
    durRest     = round(60 * ...                 % duration rest condition
                                  data_sub1.fs - 1);



%% Calc the period 
hbrSub1=data_sub1.hbr;
hbrSub2=data_sub2.hbr;

poi=[6 20];
pnoi = zeros(2,1);

  sigPart1 = [data_sub1.t, hbrSub1(:,1)];
  sigPart2 = [data_sub2.t, hbrSub2(:,1)];
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > poi(1), 1, 'first');
  pnoi(2) = find(period < poi(2), 1, 'last');
                                                                     

% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbr, 2);
coherences  = zeros(numOfChan, 2);
meanCohFP = zeros(1, length(evtTrial));                         % mean coherence in a defined spectrum for condition collaboration      
meanCohRest   = zeros(1, length(evtRest));                              % mean coherence in a defined spectrum for condition baseline
% mean coherence in a defined spectrum for condition baseline
Rsq{numOfChan} = [];
Rsq(:) = {NaN(length(period), length(data_sub1.t))};

% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hbrSub1(1, i)) && ~isnan(hbrSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    sigPart1 = [data_sub1.t, hbrSub1(:,i)];
    sigPart2 = [data_sub2.t, hbrSub2(:,i)];
    [Rsq{i}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
   
      for j=1:1:length(coi)
        Rsq{i}(period >= coi(j), j) = NaN;
      end
   
    
    % calculate mean activation in frequency band of interest
    % cooperation condition
    for j=1:1:length(evtTrial)
        try
        meanCohFP(j)  = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtTrial(j):evtTrial(j) + ...
                            durTrial)));
        catch
           meanCohFP(j) = NaN; 
        end 
    end
 

       % rest condition
    for j=1:1:length(evtRest)
        try
        meanCohRest(j)    = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtRest(j):evtRest(j) + ...
                            durRest)));
        catch
           meanCohRest(j) = NaN; 
        end 
    end
    
    fp    = nanmean(meanCohFP);                                % average mean coherences over trials
    fp_1  = nanmean(meanCohFP(1:size(meanCohFP,2)/2));
    fp_2  = nanmean(meanCohFP(size(meanCohFP,2)/2+1:end));
    rest       = nanmean(meanCohRest);


coherences(i, 1:2) = [fp,rest];
coherences1(i,1:3) = [fp_1,fp_2,rest]; 
   

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_FP_hbr_p%04d'], id, n), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'coherences','coherences1','missingtrials');
  fprintf('Data stored!\n\n');
  clear coherences
  f=1;
  catch
       
   fprintf("Repeat loop iteration");
   f=0;
  end %end for try/catch 
  end %end for while function
  end %end for permutation loop

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PD 

clear all 

srcPath = 'P:\projects\RPS\RPS\procData\hmrData\Data_PD\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_PD\';                  % processed data location
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
numOfPermutations=1000; 

for id = numOfPart
  
  % load preprocessed data
   for n=1:1:numOfPermutations
        fprintf('<strong>Permutation %d</strong>\n', n);
  f=0;
  while(f==0)
  try
  filename_sub1    = sprintf(['RPS_%02d_sub1_PD'], id);
    k=randi(numOfSources);
  filename_sub2    = sprintf(['RPS_%02d_sub2_PD'], numOfPart(k));

  fprintf('Load preprocessed data...\n');
    file_path_sub1 = strcat(srcPath, filename_sub1,'.mat');
    file_path_sub2 = strcat(srcPath, filename_sub2,'.mat');

    data_sub1=load(file_path_sub1); 
    data_sub2=load(file_path_sub2);
  %%
% -------------------------------------------------------------------------
% Check order
% -------------------------------------------------------------------------
% Column 8 indicates rest
try
evtRest     = find(data_sub1.s(:, 8) > 0);
end

% Column 10 indicates the beginning of each trial
try
evtTrial  = find(data_sub1.s(:, 10) > 0);
end

if id==18
    evtRest     = find(data_sub1.s(:, 1) > 0);
    evtTrial  = find(data_sub1.s(:, 3) > 0);
end

if size(evtTrial,1)<60
    fprintf('Trial number is different than expected 60! \n');
    missingtrials='yes';
else
    fprintf('Trial number is correct! \n');
    missingtrials='no';
end

% each trial lasts 15 seconds, rest is 60 seconds
    durTrial  = round(15 * ...                % duration task condition
                                  data_sub1.fs - 1);
    durRest     = round(60 * ...                 % duration rest condition
                                  data_sub1.fs - 1);



%% Calc the period 
hbrSub1=data_sub1.hbr;
hbrSub2=data_sub2.hbr;

poi=[6 20];
pnoi = zeros(2,1);

  sigPart1 = [data_sub1.t, hbrSub1(:,1)];
  sigPart2 = [data_sub2.t, hbrSub2(:,1)];
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > poi(1), 1, 'first');
  pnoi(2) = find(period < poi(2), 1, 'last');
                                                                     

% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbr, 2);
coherences  = zeros(numOfChan, 2);
meanCohPD = zeros(1, length(evtTrial));                         % mean coherence in a defined spectrum for condition collaboration      
meanCohRest   = zeros(1, length(evtRest));                              % mean coherence in a defined spectrum for condition baseline
% mean coherence in a defined spectrum for condition baseline
Rsq{numOfChan} = [];
Rsq(:) = {NaN(length(period), length(data_sub1.t))};

% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hbrSub1(1, i)) && ~isnan(hbrSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    sigPart1 = [data_sub1.t, hbrSub1(:,i)];
    sigPart2 = [data_sub2.t, hbrSub2(:,i)];
    [Rsq{i}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
   
      for j=1:1:length(coi)
        Rsq{i}(period >= coi(j), j) = NaN;
      end
   
    
    % calculate mean activation in frequency band of interest
    % cooperation condition
    for j=1:1:length(evtTrial)
        meanCohPD(j)  = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtTrial(j):evtTrial(j) + ...
                            durTrial)));
    end
 

       % rest condition
    for j=1:1:length(evtRest)
        try
        meanCohRest(j)    = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtRest(j):evtRest(j) + ...
                            durRest)));
        catch
           meanCohRest(j) = NaN; 
        end 
    end
    
    pd    = nanmean(meanCohPD);                                % average mean coherences over trials
    pd_1  = nanmean(meanCohPD(1:size(meanCohPD,2)/2));
    pd_2  = nanmean(meanCohPD(size(meanCohPD,2)/2+1:end));
    rest       = nanmean(meanCohRest);


coherences(i, 1:2) = [pd,rest];
coherences1(i,1:3) = [pd_1,pd_2,rest]; 
   

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_PD_hbr_p%04d'], id, n), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'coherences','coherences1','missingtrials');
  fprintf('Data stored!\n\n');
  clear coherences
  f=1;
  catch
       
   fprintf("Repeat loop iteration");
   f=0;
  end %end for try/catch 
  end %end for while function
  end %end for permutation loop

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PS 

clear all 

srcPath = 'P:\projects\RPS\RPS\procData\hmrData\Data_PS\';                        % raw data location
desPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_PS\';                  % processed data location
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
numOfPermutations=1000; 

for id = [25 26 28 29 30 31 32]
  
  % load preprocessed data
   for n=1:1:numOfPermutations
        fprintf('<strong>Permutation %d</strong>\n', n);
  f=0;
  while(f==0)
  try
  filename_sub1    = sprintf(['RPS_%02d_sub1_PS'], id);
      k=randi(numOfSources);
  filename_sub2    = sprintf(['RPS_%02d_sub2_PS'], numOfPart(k));

  fprintf('Load preprocessed data...\n');
    file_path_sub1 = strcat(srcPath, filename_sub1,'.mat');
    file_path_sub2 = strcat(srcPath, filename_sub2,'.mat');

    data_sub1=load(file_path_sub1); 
    data_sub2=load(file_path_sub2);
  %%
% -------------------------------------------------------------------------
% Check order
% -------------------------------------------------------------------------
% Column 8 indicates rest
try
evtRest     = find(data_sub1.s(:, 8) > 0);
end

% Column 10 indicates the beginning of each trial
try
evtTrial  = find(data_sub1.s(:, 10) > 0);
end

if id==18
    evtRest     = find(data_sub1.s(:, 1) > 0);
    evtTrial  = find(data_sub1.s(:, 3) > 0);
end

if size(evtTrial,1)<60
    fprintf('Trial number is different than expected 60! \n');
    missingtrials='yes';
else
    fprintf('Trial number is correct! \n');
    missingtrials='no';
end

% each trial lasts 15 seconds, rest is 60 seconds
    durTrial  = round(15 * ...                % duration task condition
                                  data_sub1.fs - 1);
    durRest     = round(60 * ...                 % duration rest condition
                                  data_sub1.fs - 1);



%% Calc the period 
hbrSub1=data_sub1.hbr;
hbrSub2=data_sub2.hbr;

poi=[6 20];
pnoi = zeros(2,1);

  sigPart1 = [data_sub1.t, hbrSub1(:,1)];
  sigPart2 = [data_sub2.t, hbrSub2(:,1)];
  [~,period,~,~,~] = wtc(sigPart1, sigPart2, 'mcc', 0); 
  pnoi(1) = find(period > poi(1), 1, 'first');
  pnoi(2) = find(period < poi(2), 1, 'last');
                                                                     

% -------------------------------------------------------------------------
% Allocate memory
% -------------------------------------------------------------------------
numOfChan = size(data_sub1.hbr, 2);
coherences  = zeros(numOfChan, 2);
meanCohPS = zeros(1, length(evtTrial));                         % mean coherence in a defined spectrum for condition collaboration      
meanCohRest   = zeros(1, length(evtRest));                              % mean coherence in a defined spectrum for condition baseline
% mean coherence in a defined spectrum for condition baseline
Rsq{numOfChan} = [];
Rsq(:) = {NaN(length(period), length(data_sub1.t))};

% -------------------------------------------------------------------------
% Calculate Coherence increase between conditions for every channel of the 
% dyad
% -------------------------------------------------------------------------
fprintf('<strong>Estimation of the wavelet transform coherence for all channels...</strong>\n');
for i=1:1:numOfChan
  if ~isnan(hbrSub1(1, i)) && ~isnan(hbrSub2(1, i))                         % check if this channel was not rejected in both subjects during preprocessing
    sigPart1 = [data_sub1.t, hbrSub1(:,i)];
    sigPart2 = [data_sub2.t, hbrSub2(:,i)];
    [Rsq{i}, ~, ~, coi, ~] = wtc(sigPart1, sigPart2, 'mcc', 0);                % r square - measure for coherence
  
   
      for j=1:1:length(coi)
        Rsq{i}(period >= coi(j), j) = NaN;
      end
   
    
    % calculate mean activation in frequency band of interest
    % cooperation condition
    for j=1:1:length(evtTrial)
        meanCohPS(j)  = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtTrial(j):evtTrial(j) + ...
                            durTrial)));
    end
 

       % rest condition
    for j=1:1:length(evtRest)
        try
        meanCohRest(j)    = nanmean(nanmean(Rsq{i}(pnoi(1):pnoi(2), ...
                            evtRest(j):evtRest(j) + ...
                            durRest)));
        catch
           meanCohRest(j) = NaN; 
        end 
    end
    
    ps    = nanmean(meanCohPS);                                % average mean coherences over trials
    ps_1  = nanmean(meanCohPS(1:size(meanCohPS,2)/2));
    ps_2  = nanmean(meanCohPS(size(meanCohPS,2)/2+1:end));
    rest       = nanmean(meanCohRest);


coherences(i, 1:2) = [ps,rest];
coherences1(i,1:3) = [ps_1,ps_2,rest]; 
   

end

  
end
% save preprocessed data
  desFolder   = strcat(desPath);
  
  
  file_path = strcat(desFolder, sprintf(['RPS_%02d_PS_hbr_p%04d'], id, n), ...
                     '.mat');

  fprintf('The wtc data of dyad will be saved in '); 
  fprintf('%s ...\n', file_path);
  save(file_path, 'coherences','coherences1','missingtrials');
  fprintf('Data stored!\n\n');
  clear coherences
  f=1;
  catch
       
   fprintf("Repeat loop iteration");
   f=0;
  end %end for try/catch 
  end %end for while function
  end %end for permutation loop

end