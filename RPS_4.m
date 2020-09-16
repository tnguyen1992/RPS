%%%%%%%%% Extract WTC data %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
clear all 

srcPath = 'P:\projects\RPS\RPS\procData\wtcData\Data_C\';                  % source data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_C.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_C.mat'));
  end
  %%
  for i=numOfPart
      
    
    filename    = sprintf(['RPS_%02d_C'], i);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
    
for cond=1
  % dyad ID
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),2)  = [cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond];
  % ch
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),3)  = [1:16];
  %roi
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),5)  =  coherences(1:16,cond);
end
  for rest=1
  rest_coherences(:,rest,i)=coherences(:,2);
  end
  
  end
  
  
 srcPath = 'P:\projects\RPS\RPS\procData\wtcData\Data_FP\';                  % source data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_FP.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_FP.mat'));
  end
  %%
  for i=numOfPart
      
    
    filename    = sprintf(['RPS_%02d_FP'], i);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
    
for cond=3
  % dyad ID
  % dyad ID
  wtc((length(coherences)*(cond-1)-15+((i-1)*80)):...
      (length(coherences)*(cond-1)+((i-1)*80)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  wtc((length(coherences)*(cond-1)-15+((i-1)*80)):...
      (length(coherences)*(cond-1)+((i-1)*80)),2)  = [cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1];
   % ch
  wtc((length(coherences)*(cond-1)-15+((i-1)*80)):...
      (length(coherences)*(cond-1)+((i-1)*80)),3)  = [1:16];
  %roi
  wtc((length(coherences)*(cond-1)-15+((i-1)*80)):...
      (length(coherences)*(cond-1)+((i-1)*80)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  wtc((length(coherences)*(cond-1)-15+((i-1)*80)):...
      (length(coherences)*(cond-1)+((i-1)*80)),5)  =  coherences(1:16,cond-2);
  end
   for rest=2
  rest_coherences(:,rest,i)=coherences(:,2);
  end 
  end 
  
  srcPath = 'P:\projects\RPS\RPS\procData\wtcData\Data_PD\';                  % source data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_PD.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_PD.mat'));
  end
  %%
  for i=numOfPart
      
    
    filename    = sprintf(['RPS_%02d_PD'], i);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
    
for cond=5
  % dyad ID
  % dyad ID
  wtc((length(coherences)*(cond-2)-15+((i-1)*80)):...
      (length(coherences)*(cond-2)+((i-1)*80)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  wtc((length(coherences)*(cond-2)-15+((i-1)*80)):...
      (length(coherences)*(cond-2)+((i-1)*80)),2)  = [cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2 cond-2];
  % ch
  wtc((length(coherences)*(cond-2)-15+((i-1)*80)):...
      (length(coherences)*(cond-2)+((i-1)*80)),3)  = [1:16];
  %roi
  wtc((length(coherences)*(cond-2)-15+((i-1)*80)):...
      (length(coherences)*(cond-2)+((i-1)*80)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  wtc((length(coherences)*(cond-2)-15+((i-1)*80)):...
      (length(coherences)*(cond-2)+((i-1)*80)),5)  =  coherences(1:16,cond-4);
  end
  for rest=3
  rest_coherences(:,rest,i)=coherences(:,2);
  end
  end
  
  srcPath = 'P:\projects\RPS\RPS\procData\wtcData\Data_PS\';                  % source data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*_PS.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_PS.mat'));
  end
  %%
  for i=numOfPart
      
    
    filename    = sprintf(['RPS_%02d_PS'], i);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
    
for cond=7
  % dyad ID
  % dyad ID
  wtc((length(coherences)*(cond-3)-15+((i-1)*80)):...
      (length(coherences)*(cond-3)+((i-1)*80)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  wtc((length(coherences)*(cond-3)-15+((i-1)*80)):...
      (length(coherences)*(cond-3)+((i-1)*80)),2)  = [cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3 cond-3];
   % ch
  wtc((length(coherences)*(cond-3)-15+((i-1)*80)):...
      (length(coherences)*(cond-3)+((i-1)*80)),3)  = [1:16];
  %roi
  wtc((length(coherences)*(cond-3)-15+((i-1)*80)):...
      (length(coherences)*(cond-3)+((i-1)*80)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  wtc((length(coherences)*(cond-3)-15+((i-1)*80)):...
      (length(coherences)*(cond-3)+((i-1)*80)),5)  =  coherences(1:16,cond-6);
  end
  for rest=4
  rest_coherences(:,rest,i)=coherences(:,2);
  end  
  end
 
  
 %% average all coherences 
rest_coher=mean(rest_coherences,2);
  
  for i=1:32
      
  cond=5
  % dyad ID
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),2)  = [cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond];
  % ch
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),3)  = [1:16];
  %roi
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  wtc((length(coherences)*cond-15+((i-1)*80)):...
      (length(coherences)*cond+((i-1)*80)),5)  =  rest_coher(1:16,1,i);
end
  
  
%% %%%%%%% Extract GLM data %%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%
% clear all 
% 
% srcPath = '\\fs.univie.ac.at\homedirs\nguyenq22\Documents\MATLAB\scripts\wieki_fnirs\procData\glmData\';                  % source data location
% %% Scan for all subjects
%   sourceList    = dir([srcPath, '*.mat']);
%   sourceList    = struct2cell(sourceList);
%   sourceList    = sourceList(1,:);
%   numOfSources  = length(sourceList);
%   PartList       = zeros(1, numOfSources);
%   PartList  = erase(sourceList, '.mat');
% 
%   for i=1:1:numOfSources
%       
%     srcFolder   = strcat(srcPath);
%     filename    = PartList{1, i};
% 
%     fprintf('Load preprocessed data...\n');
%     file_path = strcat(srcFolder, filename,'.mat');
%     load(file_path);
%     
%     for cond=1:4 
%   % dyad ID
%   glm((length(beta)*cond-18+((i-1)*79)):...
%       (length(beta)*cond+((i-1)*79)),1)  =  [i i i i i i i i i i i i i i i i i i i];
% %   wtc((length(coherences)*cond-21+((i-1)*88)):...
% %       (length(coherences)*cond+((i-1)*88)),2)  =  [filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename filename];
%   % cond
%   glm((length(beta)*cond-18+((i-1)*79)):...
%       (length(beta)*cond+((i-1)*79)),2)  = [cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond cond];
%   % ch
%   glm((length(beta)*cond-18+((i-1)*79)):...
%       (length(beta)*cond+((i-1)*79)),3)  = [1:19];
%   % roi (1 tpj) (0 NaN) (2  dlpfc) (3 superior) (4 frontopolar) (5
%   % inferior)
%   glm((length(beta)*cond-18+((i-1)*79)):...
%       (length(beta)*cond+((i-1)*79)),4)  = [1 2 2 1 1 5 2 3 1 5 3 1 1 5 4 3 1 4 4];
%   % wtc
%   glm((length(beta)*cond-18+((i-1)*79)):...
%       (length(beta)*cond+((i-1)*79)),5)  =  beta(1:19,cond);
%   end
%     
%   end
  