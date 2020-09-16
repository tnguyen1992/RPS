srcPath = 'P:\projects\RPS\RPS\procData\gcData\Data_FP\';                  % source data location
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
    
    sign(:,1) = granger(:,1)-granger(:,2);
    sign(:,2) = granger(:,3)-granger(:,4);
    
    for ch=1:16
    if (sign(ch,1) > 0) && (sign(ch,2) > 0)
        sign(ch,3)=3;
    elseif (sign(ch,1) > 0) && (sign(ch,2) < 0)
        sign(ch,3)=1;
    elseif (sign(ch,1) < 0) && (sign(ch,2) > 0)
        sign(ch,3)=2;
    else
        sign(ch,3)=0;
    end
    end
    
    for cond=3
  % dyad ID
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),2)  = [cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1];
   % ch
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),3)  = [1:16];
  %roi
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),5)  =  sign(:,3);
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
    
    for cond=4
 % dyad ID
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),2)  = [cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1];
   % ch
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),3)  = [1:16];
  %roi
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),5)  =  sign(:,3);
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
    
for cond=5
    % dyad ID
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),1)  =  [i i i i i i i i i i i i i i i i];
  % cond
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),2)  = [cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1 cond-1];
   % ch
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),3)  = [1:16];
  %roi
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),4)  = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4];
  % wtc
  gc_sign((length(sign)*(cond-1)-15+((i-1)*48)):...
      (length(sign)*(cond-1)+((i-1)*48)),5)  =  sign(:,3);
    end
  end
    
    dlmwrite("P:\projects\RPS\RPS\gc_df.csv",gc_sign);
  
    
    