%%%%%%%%% Extract WTC data %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
clear all 

srcPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_C\';                  % source data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_*.mat'));
  end
  %%

 
  for i=numOfPart
      for rpa_id=1:1:1000

    
    filename    = sprintf(['RPS_%02d_C_p%04d'], i, rpa_id);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
      rpa_c(1:16,1:2,rpa_id)=coherences;
      
      
      end
rpa=nanmean(rpa_c,3);

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
      (length(coherences)*cond+((i-1)*80)),5)  =  rpa(1:16,cond);
end
  for rest=1
  rest_coherences(:,rest,i)=rpa(:,2);
  end
  
  end
  
  clear rpa
 srcPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_FP\';                  % source data location
%% Scan for all subjects
  sourceList    = dir([srcPath, 'RPS_*.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_*.mat'));
  end
  %%
  for i=numOfPart
      
       for rpa_id=1:1:1000

    filename    = sprintf(['RPS_%02d_FP_p%04d'], i, rpa_id);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
    rpa_c(1:16,1:2,rpa_id)=coherences;
      
      
      end
rpa=nanmean(rpa_c,3);
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
      (length(coherences)*(cond-1)+((i-1)*80)),5)  =  rpa(1:16,cond-2);
  end
   for rest=2
  rest_coherences(:,rest,i)=rpa(:,2);
  end 
  end 
  
  clear rpa
  srcPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_PD\';                  % source data location
%% Scan for all subjects
   sourceList    = dir([srcPath, 'RPS_*.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_*.mat'));
  end
  %%
  for i=numOfPart
             for rpa_id=1:1:1000
 
    
    filename    = sprintf(['RPS_%02d_PD_p%04d'], i, rpa_id);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
     rpa_c(1:16,1:2,rpa_id)=coherences;
      
      
      end
rpa=nanmean(rpa_c,3);
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
      (length(coherences)*(cond-2)+((i-1)*80)),5)  =  rpa(1:16,cond-4);
  end
  for rest=3
  rest_coherences(:,rest,i)=rpa(:,2);
  end
  end
  
  clear rpa
  srcPath = 'P:\projects\RPS\RPS\procData\rpaData\Data_PS\';                  % source data location
%% Scan for all subjects
   sourceList    = dir([srcPath, 'RPS_*.mat']);
  sourceList    = struct2cell(sourceList);
  sourceList    = sourceList(1,:);
  numOfSources  = length(sourceList);
  PartList       = zeros(1, numOfSources);
  for ii=1:1:numOfSources
    numOfPart(ii)  = sscanf(sourceList{ii}, ...
                    strcat('RPS_%d_*.mat'));
  end
  %%
  for i=numOfPart
                   for rpa_id=1:1:1000

    
    filename    = sprintf(['RPS_%02d_PS_p%04d'], i, rpa_id);

    fprintf('Load preprocessed data...\n');
    file_path = strcat(srcPath, filename,'.mat');
    load(file_path);
    rpa_c(1:16,1:2,rpa_id)=coherences;
      
      
      end
rpa=nanmean(rpa_c,3);
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
      (length(coherences)*(cond-3)+((i-1)*80)),5)  =  rpa(1:16,cond-6);
  end
  for rest=4
  rest_coherences(:,rest,i)=rpa(:,2);
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
  
  
