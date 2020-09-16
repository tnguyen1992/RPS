srcPath = 'P:\projects\RPS\RPS\rawData\';
sourceList    = dir([srcPath, 'RPS_*']);
sourceList    = struct2cell(sourceList);
sourceList    = sourceList(1,:);
% wie viele Probanden gibt es?
numOfSources  = length(sourceList);
numOfPart       = zeros(1, numOfSources);

  for i=1:1:numOfSources
    numOfPart(i)  = sscanf(sourceList{i}, ['RPS_%02d']);
  end
  
  for i = numOfPart
  srcFolder   = strcat(srcPath, sprintf(['RPS_%02d/'], i));
  
  Sub1SrcDir = strcat(srcFolder, sprintf(['Subject1/']));
  Sub2SrcDir = strcat(srcFolder, sprintf(['Subject2/']));
  
% -------------------------------------------------------------------------
% Wo sind die Daten im jeweiligen Ordner? 
% -------------------------------------------------------------------------

  Sub1_NIRSFile = strcat(Sub1SrcDir, sprintf(['RPS_%02d_FP.mat'], i));
  load(Sub1_NIRSFile);
  hbo=Y.hbo;
for j = 1:1:size(hbo, 2)
  subplot(4,4,j);
  fs=7.8125;
  t = 0:1/fs:length(hbo)/fs - 1/fs;
    t = t';                                                                     %#ok<NASGU>

  sig = [t, hbo(:,j)];
  sigma2=var(sig(:,2));                                                     % estimate signal variance
  
  [wave,period,~,coi,~] = wt(sig);                                          % compute wavelet power spectrum
  power = (abs(wave)).^2 ;
  
  for k=1:1:length(coi)
    wave(period >= coi(k), k) = NaN;                                        % set values below cone of interest to NAN
  end

     h = imagesc(t, log2(period), log2(abs(power/sigma2)));
    colorbar;
    Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
    set(gca,'YLim',log2([min(period),max(period)]), ...
          'YDir','reverse', 'layer','top', ...
          'YTick',log2(Yticks(:)), ...
          'YTickLabel',num2str(Yticks'), ...
          'layer','top')
    title(sprintf('Channel %d', j));
    ylabel('Period in seconds');
    xlabel('Time in seconds');
    set(h, 'AlphaData', ~isnan(wave));

    colormap jet;
    end
    set(gcf,'units','normalized','outerposition',[0 0 1 1])                     % maximize figure

    badChannels = RPS_channelCheckbox();
    badChannels_collect{i,1}=badChannels;
    close(gcf); 

  Sub2_NIRSFile = strcat(Sub2SrcDir, sprintf(['RPS_%02d_FP.mat'], i));
  
  load(Sub2_NIRSFile);
  hbo=Y.hbo;
for j = 1:1:size(hbo, 2)
  subplot(4,4,j);
  fs=7.8125;
  t = 0:1/fs:length(hbo)/fs - 1/fs;
    t = t';  
  sig = [t, hbo(:,j)];
  sigma2=var(sig(:,2));                                                     % estimate signal variance
  
  [wave,period,~,coi,~] = wt(sig);                                          % compute wavelet power spectrum
  power = (abs(wave)).^2 ;
  
  for k=1:1:length(coi)
    wave(period >= coi(k), k) = NaN;                                        % set values below cone of interest to NAN
  end

     h = imagesc(t, log2(period), log2(abs(power/sigma2)));
    colorbar;
    Yticks = 2.^(fix(log2(min(period))):fix(log2(max(period))));
    set(gca,'YLim',log2([min(period),max(period)]), ...
          'YDir','reverse', 'layer','top', ...
          'YTick',log2(Yticks(:)), ...
          'YTickLabel',num2str(Yticks'), ...
          'layer','top')
    title(sprintf('Channel %d', i));
    ylabel('Period in seconds');
    xlabel('Time in seconds');
    set(h, 'AlphaData', ~isnan(wave));

    colormap jet;
    end
    set(gcf,'units','normalized','outerposition',[0 0 1 1])                     % maximize figure

    badChannels = RPS_channelCheckbox();
    badChannels_collect{i,2}=badChannels;
    close(gcf); 
  
  end
  