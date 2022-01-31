
% The following code has been used to analyze the data presented in the 
% paper:
%'Low-frequency subthalamic neural oscillations are involved in explicit 
% and implicit facial emotional processing - a local field potential study'
% (Duprez et al., XXXX)
% Time-frequency analyses are done using complex morlet wavelet convolution
% based on published equations (Cohen 2014 - Analyzing neural time series 
% data).
% This code has been adapted from Mike X Cohen's one.
% Although at many places efficiency of this code can be improved, 
% it does correct computations.


%% LFPemo analysis

% loaded files have are an EEG structure such as in eeglab
matfiles = dir(fullfile(['*.mat']));
numfiles = length(matfiles);

%% ERP
ERP = zeros(numfiles, 3, 3000);
for subno = 1:numfiles
    
    load(matfiles(subno).name) % load file
    LFPemo.data = double(LFPemo.data); % convert to double precision
    
    for chani = 1:3
    ERP(subno, chani, :) = squeeze(mean(LFPemo.data(chani,:,:), 3));
    end
end

contact = {'0-1', '1-2', '2-3'};
figure; 
for chani = 1:3
    subplot(3, 1, chani)
    plot(LFPemo.times(501:3000),squeeze(mean(ERP(:, chani, 501:3000), 1)))
    set(gca, 'ylim', [-0.6, 0.6]);
    title(['ERP - Contact pair ' num2str(contact{chani})])
end

%% Time-frequency analyses
% Soft-coded convolution parameters

min_frex  =  1;
max_frex  = 40;
num_frex  = 50;

% set range for number of wavelet cycles
nCycRange = [4 10];

times2save = -500:15:2000; % downsample to 100 Hz
basetime   = [-600 -300]; % baseline period

% Load data
% initialize output time-frequency data
% 4 for the number of conditions, 3 number of channels
tf = zeros(numfiles, 4, 3, num_frex,size(times2save,2)); 


% initialize output itpc data
% 4 for the number of conditions, 3 number of channels
itpcz = zeros(numfiles, 4,3,num_frex,size(times2save,2)); 


for subno = 1:numfiles
    
    load(matfiles(subno).name)
    LFPemo.data = double(LFPemo.data); % convert to double precision
    
    % Here, remove phase-locked activity (ERP subtraction)
    
    LFPemo.data = LFPemo.data - squeeze(mean(LFPemo.data, 3));
    
    % With data reflection to avoid edge artifacts

    LFPemo.data_reflect = [LFPemo.data(:,end:-1:1,:) LFPemo.data(:,:,:) ...
        LFPemo.data(:,end:-1:1,:) ];
    % Set times2save for temporal downsampling at 67 Hz
    
    times2saveidx = dsearchn(LFPemo.times',times2save');
    
    frex = logspace(log10(min_frex),log10(max_frex),num_frex);
    
    % wavelet parameters
    s = logspace(log10(nCycRange(1)),log10(nCycRange(2)),...
        num_frex)./(2*pi.*frex);
    t = -2:1/LFPemo.srate:2;
    halfwave = floor((length(t)-1)/2);
    
    % Data convolution parameters
    nData = size(LFPemo.data_reflect, 3)*size(LFPemo.data_reflect, 2);
    nWave = length(t);
    nConv = nData+nWave-1;
    
    % Data wavelets
    cmwX = zeros(num_frex,nConv);
    for fi=1:num_frex
        cmw = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2))...
            ,nConv);
        % amplitude-normalize in the frequency domain
        cmwX(fi,:) = cmw ./ max(cmw); 
    end
    
    % indices for the baseline
    baseidx = dsearchn(LFPemo.times',basetime');
    
    % loop through channels
    for chani = 1:size(LFPemo.data_reflect,1)
        dataX = fft(reshape(LFPemo.data_reflect(chani,:,:),1,nData)...
            ,nConv); % concatenates all the trials of the first channel
        
        for fi=1:num_frex
            as = ifft(dataX.*cmwX(fi,:) );
            as = as(halfwave+1:end-halfwave);
            as = reshape(as,1, size(LFPemo.data_reflect, 2), ...
                size(LFPemo.data_reflect, 3));
            as = as(:, size(LFPemo.data, 2)+1:size(LFPemo.data, 2)*2,:);
            
            % condition-average baseline
            basepow = mean(squeeze(mean(abs(as(:,baseidx(1):baseidx(2),...
                :)).^2,2)),1);
            baseitpc = mean(abs(mean(squeeze(exp(1i*angle(as(:,...
                baseidx(1):baseidx(2),:)))),2)),1);
            
            
            for condi = 1:4
                tf(subno,condi,chani,fi,:) = 10*log10(mean(abs(as(:,...
                    times2saveidx,LFPemo.condition  == condi)).^2,3)...
                    ./ basepow ); 
                itpcz(subno,condi,chani,fi,:) = size(as(:,:,...
                    LFPemo.condition== condi),3)*...
                    (10*log10(abs(mean(exp(1i*angle(as(:,...
                    times2saveidx,LFPemo.condition == condi))),3))...
                    ./ baseitpc).^2);
                
            end % end condition loop
            
        end % end frequencies loop
        
    end% end channel loop
    
    disp(['subject ', num2str(subno), ' is done'])
end% subject loop



%% Plot

% % Overall power
%
%

figure(24), clf

climdb  = [-1 1];
contourf(times2save,frex,squeeze(mean(mean(squeeze(tf(:,:,1,:,:))...
    ,2) ,1)), 20,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Averaged time-frequency power')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')
cbh.FontSize = 28
% ylabel(colorbar, 'dB  from baseline', 'FontSize', 18)
ax = gca;
ax.FontSize = 20;
x0=10;
y0=10;
width=800;
height=500;

%
% Power by condition
figure(2), clf
condiname = {'Gender-Neutral'; 'Gender-Fear'; 'Emotion-Neutral';...
    'Emotion-Fear'; 'Congruent 1€'; 'Incongruent 1€'};
climdb  = [-1 1];
for condi=1:4
    subplot(2,2,condi)
    contourf(times2save,frex,squeeze(mean(squeeze(tf(:,condi,1,:,:))...
        ,1)),20,'linecolor','none')
    set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex)...
        ,log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
        log10(max_frex),6)*10)/10,'xlim', [-300 2000])
    title(condiname(condi))
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    colorbar
    ylabel(colorbar, 'dB  from baseline')
    %colormap jet
end
%

%% Permutation testing to isolate significant clusters of tf power 
%  differences from baseline
load(matfiles(1).name) % to get the time information from one dataset

% take only first bipolar montage (channel 1)
tfavg = squeeze(mean(tf(:,:,1,:,:), 2)); 

nTimepoints = numel(times2save);
baseidx = dsearchn(times2save',[-500 -200]');

voxel_pval   = 0.05;
cluster_pval = 0.05;
n_permutes = 1000;

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_frex);
permuted_vals    = zeros(n_permutes,size(tfavg, 1), ...
    num_frex,numel(times2save));
max_clust_info   = zeros(n_permutes,1);


for permi=1:n_permutes
    
    for subi = 1:size(tfavg, 1)
        cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1);
        permuted_vals(permi,subi,:,:) = tfavg(subi,:,[cutpoint:end ...
            1:cutpoint-1]);
    end
    disp(['permutation ', num2str(permi), ' is done']);
end
realmean = squeeze(mean(tfavg, 1));
perm_mean = squeeze(mean(permuted_vals, 2));
zmap = (realmean-squeeze(mean(perm_mean,1)))./squeeze(std(perm_mean,1));

threshmean = realmean;
threshmean(abs(zmap)<norminv(1-voxel_pval))=0;



figure(28), clf
subplot(221)
climdb  = [-1 1];
contourf(times2save,frex,realmean, 60,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Averaged time-frequency power')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')
colormap jet


subplot(222)
climdb  = [-1 1];
contourf(times2save,frex,zmap, 60,'linecolor','none')
set(gca,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')
colormap jet

subplot(223)
contourf(times2save,frex,threshmean, 60,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Uncorrected power map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')
colormap jet



% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    % for cluster correction, apply uncorrected threshold and 
    % get maximum cluster sizes
    fakecorrsz = squeeze((perm_mean(permi,:,:)-mean(perm_mean,1))...
        ./ std(perm_mean,[],1) );
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(fakecorrsz);
    max_clust_info(permi) = max([ 0 cellfun(@numel,...
        clustinfo.PixelIdxList)]); % zero accounts for empty maps
    % using cellfun here eliminates the need for a slower loop over cells
end

% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end


% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<clust_threshold
        zmap(islands.PixelIdxList{i})=0;
    end
end




zmapthresh=logical(zmapthresh);




figure
climdb  = [-1, 1];
contourf(times2save,frex,realmean, 20,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Averaged time-frequency power')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')
%colormap jet
hold on
contour(times2save,frex,zmapthresh, 1,'linecolor','k')
ylabel('Frequency (Hz)')
xlabel('Time (ms)')
cbh.Label.String = 'dB  from baseline'
cbh.FontSize = 28
% ylabel(colorbar, 'dB  from baseline', 'FontSize', 18)
ax = gca;
ax.FontSize = 20;
x0=10;
y0=10;
width=800;
height=500;
print('figure', '-dpng', '-r1000')


%% Extract power from the tf windows of interest

time = times2save;

%% Power

%% Delta
% Definition of TF windows

time_windows = [ 490 760 ]; % 3 time windows
freq_windows = [ frex(11) frex(18) ]; % 3 frequency windows

% find indices corresponding to time and frequency windows

timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));

for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(times2save-time_windows(i,j)));
    end
end

%for the frequency windows
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(frex-freq_windows(i,j)));
    end
end

for condi=1:size(tf,2) %
    
    % pointer to stats file
    statsfilename=(['statistics_file_TF_delta_reflected_rmtl', num2str(condi), '.txt']);
    
    fid=fopen(statsfilename,'w');
    
    for subno=0:size(tf,1)
        
        
        % Write out column variable name or subject number.
        if subno==0
            fprintf(fid,'subnum\t');
        else
            fprintf(fid,'%g\t',subno);
        end
        
        for ti=1:size(time_windows,1)
            
            for fi=1:size(freq_windows,1)
                
                % Write out the column variable name of data.
                if subno==0
                    
                    fprintf(fid,[ num2str(time_windows(ti,1)) ...
                        num2str(time_windows(ti,2)) '_f ' ...
                        num2str(freq_windows(fi,1)) ...
                        num2str(freq_windows(fi,2)) '\t' ]);
                    
                else
                    
                    % get data from large window
                    % Note: This is condition-specific; you could make 
                    % it condition-average by taking the mean across 
                    % the 2nd dimension.
                    tempdat = squeeze(tf(subno,condi,1,...
                        freqidx(fi,1):freqidx(fi,2),...
                        timeidx(ti,1):timeidx(ti,2)));
                    
                    % find max point (must first vectorize otherwise you 
                    % get max for each line)
                    [junk,maxpoint]=max(tempdat(:));
                    
                    % convert that back to 2d coordinates
                    [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                    
                    % get the indices from the larger matrix, not the 
                    % indices of the smaller TF window.
                    maxF = maxF+freqidx(1)-1;
                    maxT = maxT+timeidx(1)-1;
                    
                    % Now write out data and max time/frequency to file.
                    % 6 so as to have a proportional number of frequency 
                    % (3) and timepoints in the subject TF window
                    fprintf(fid,'%g\t',mean(mean(tempdat))); 
                end
                
            end
            
        end
        
        % CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
        fprintf(fid,'\n'); 
        if subno > 0
            maxfreqdelta1(subno)=frex(maxF);
        else
            continue
        end
    end % end subject loop
    
    fclose(fid);
end % end condition loop



%% Alpha
% Definition of TF windows

time_windows = [ 1435 1780 ]; % 3 time windows
freq_windows = [ frex(29) frex(35) ]; % 3 frequency windows

% find indices corresponding to time and frequency windows

timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));

for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(times2save-time_windows(i,j)));
    end
end

%for the frequency windows
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(frex-freq_windows(i,j)));
    end
end

for condi=1:size(tf,2) %
    
    % pointer to stats file
    statsfilename=(['statistics_file_TF_alpha_reflected_rmtl_', num2str(condi), '.txt']);
    
    fid=fopen(statsfilename,'w');
    
    for subno=0:size(tf,1)
        
        
        % Write out column variable name or subject number.
        if subno==0
            fprintf(fid,'subnum\t');
        else
            fprintf(fid,'%g\t',subno);
        end
        
        for ti=1:size(time_windows,1)
            
            for fi=1:size(freq_windows,1)
                
                % Write out the column variable name of data.
                if subno==0
                    
                    fprintf(fid,[ num2str(time_windows(ti,1)) ...
                        num2str(time_windows(ti,2)) '_f ' ...
                        num2str(freq_windows(fi,1)) ...
                        num2str(freq_windows(fi,2)) '\t' ]);
                    
                else
                    
                    % get data from large window
                    % Note: This is condition-specific; you could make 
                    % it condition-average by taking the mean across the 
                    % 2nd dimension.
                    tempdat = squeeze(tf(subno,condi,1,freqidx(fi,1):...
                        freqidx(fi,2),timeidx(ti,1):timeidx(ti,2)));
                    
                    % find max point (must first vectorize otherwise 
                    % you get max for each line)
                    [junk,maxpoint]=max(tempdat(:));
                    
                    % convert that back to 2d coordinates
                    [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                    
                    % get the indices from the larger matrix, not the 
                    %indices of the smaller TF window.
                    maxF = maxF+freqidx(1)-1;
                    maxT = maxT+timeidx(1)-1;
                    
                    % Now write out data and max time/frequency to file.
                    % 6 so as to have a proportional number of frequency 
                    % (3) and timepoints in the subject TF window
                    fprintf(fid,'%g\t',mean(mean(tempdat))); 
                end
            end
        end
        
        % CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
        fprintf(fid,'\n'); 
        if subno > 0
            maxfreqdelta2(subno)=frex(maxF);
        else
            continue
        end
    end % end subject loop
    
    fclose(fid);
end % end condition loop



%% Beta
% Definition of TF windows

time_windows = [ 1390 1810 ]; % 3 time windows
freq_windows = [ frex(37) frex(46) ]; % 3 frequency windows

% find indices corresponding to time and frequency windows

timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));

for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(times2save-time_windows(i,j)));
    end
end

%for the frequency windows
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(frex-freq_windows(i,j)));
    end
end

for condi=1:size(tf,2) %
    
    % pointer to stats file
    statsfilename=(['statistics_file_TF_beta_reflected_rmtl_', num2str(condi), '.txt']);
    
    fid=fopen(statsfilename,'w');
    
    for subno=0:size(tf,1)
        
        
        % Write out column variable name or subject number.
        if subno==0
            fprintf(fid,'subnum\t');
        else
            fprintf(fid,'%g\t',subno);
        end
        
        for ti=1:size(time_windows,1)
            
            for fi=1:size(freq_windows,1)
                
                % Write out the column variable name of data.
                if subno==0
                    
                    fprintf(fid,[ num2str(time_windows(ti,1))...
                        num2str(time_windows(ti,2)) '_f '...
                        num2str(freq_windows(fi,1)) ...
                        num2str(freq_windows(fi,2)) '\t' ]);
                    
                else
                    
                    % get data from large window
                    % Note: This is condition-specific; 
                    % you could make it condition-average by 
                    % taking the mean across the 2nd dimension.
                    tempdat = squeeze(tf(subno,condi,1,freqidx(fi,1):...
                        freqidx(fi,2),timeidx(ti,1):timeidx(ti,2)));
                    
                    % find max point (must first vectorize otherwise 
                    % you get max for each line)
                    [junk,maxpoint]=min(tempdat(:));
                    
                    % convert that back to 2d coordinates
                    [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                    
                    % get the indices from the larger matrix, 
                    % not the indices of the smaller TF window.
                    maxF = maxF+freqidx(1)-1;
                    maxT = maxT+timeidx(1)-1;
                    
                    % Now write out data and max time/frequency to file.
                    % 6 so as to have a proportional number of frequency 
                    % (3) and timepoints in the subject TF window
                    fprintf(fid,'%g\t',mean(mean(tempdat))); 
                end
                
            end
            
        end
        
        % CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
        fprintf(fid,'\n');
        if subno > 0
            maxfreqbeta(subno)=frex(maxF);
        else
            continue
        end
    end % end subject loop
    
    fclose(fid);
end % end condition loop

maxfreqdelta1=maxfreqdelta1';
maxfreqdelta2=maxfreqdelta2';
maxfreqbeta=maxfreqbeta';
save('maxFreq_delta1.txt', 'maxfreqdelta1', '-ascii', '-double', '-tabs')
save('maxFreq_delta2.txt', 'maxfreqdelta2', '-ascii', '-double', '-tabs')
save('maxFreq_beta.txt', 'maxfreqbeta', '-ascii', '-double', '-tabs')


%% Same thing for ITPCZ

%% Plot

figure(23), clf

climdb  = [0 200];
contourf(times2save,frex,squeeze(mean(mean(squeeze(itpcz(:,:,1,:,:))...
    ,2) ,1)),20,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Averaged ITPCz')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'ITPCz')
cbh.FontSize = 28
% ylabel(colorbar, 'dB  from baseline', 'FontSize', 18)
ax = gca;
ax.FontSize = 20;
x0=10;
y0=10;
width=800;
height=500;

%
%
% ITPCz by condition
figure(2), clf
condiname = {'Gender-Neutral'; 'Gender-Fear'; 'Emotion-Neutral';...
    'Emotion-Fear'; 'Congruent 1€'; 'Incongruent 1€'};
climdb  = [0 250];
for condi=1:4
    subplot(2,2,condi)
    contourf(times2save,frex,squeeze(mean(squeeze(itpcz(:,condi,1,:,:))...
        ,1)),60,'linecolor','none')
    set(gca,'clim',climdb, 'yscale','log','ytick',logspace(log10(min_frex),...
        log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
        log10(max_frex),6)*10)/10,'xlim', [-300 2000])
    title(condiname(condi))
    xlabel('Time (s)'), ylabel('Frequency (Hz)')
    colorbar
    ylabel(colorbar, 'ITPCz')
    % colormap jet
end
%
%

%% Permutation testing to isolate significant clusters 
%  of tf power differences from baseline
load(matfiles(1).name) % to get the time information from one dataset

% take only first bipolar montage (channel 1)
tfavg = squeeze(mean(itpcz(:,:,1,:,:), 2)); 

nTimepoints = numel(times2save);
baseidx = dsearchn(times2save',[-500 -200]');

voxel_pval   = 0.05;
cluster_pval = 0.05;
n_permutes = 1000;

% initialize null hypothesis matrices
permuted_maxvals = zeros(n_permutes,2,num_frex);
permuted_vals    = zeros(n_permutes,size(tfavg, 1),...
    num_frex,numel(times2save));
max_clust_info   = zeros(n_permutes,1);


for permi=1:n_permutes
    
    for subi = 1:size(tfavg, 1)
        cutpoint = randsample(2:nTimepoints-diff(baseidx)-2,1);
        permuted_vals(permi,subi,:,:) = tfavg(subi,:,...
            [cutpoint:end 1:cutpoint-1]);
    end
    disp(['permutation ', num2str(permi), ' is done']);
end
realmean = squeeze(mean(tfavg, 1));
perm_mean = squeeze(mean(permuted_vals, 2));
zmap = (realmean-squeeze(mean(perm_mean,1)))./squeeze(std(perm_mean,1));

threshmean = realmean;
threshmean(abs(zmap)<norminv(1-voxel_pval))=0;


figure(28), clf
subplot(221)
contourf(times2save,frex,realmean, 60,'linecolor','none')
set(gca,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Averaged time-frequency power')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')


subplot(222)
contourf(times2save,frex,zmap, 60,'linecolor','none')
set(gca,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('unthresholded Z map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')

subplot(223)
contourf(times2save,frex,threshmean, 60,'linecolor','none')
set(gca,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Uncorrected power map')
xlabel('Time (ms)'), ylabel('Frequency (Hz)')
colorbar
ylabel(colorbar, 'dB  from baseline')



% this time, the cluster correction will be done on the permuted data, thus
% making no assumptions about parameters for p-values
for permi = 1:n_permutes
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    fakecorrsz = squeeze((perm_mean(permi,:,:)-mean(perm_mean,1))...
        ./ std(perm_mean,[],1) );
    fakecorrsz(abs(fakecorrsz)<norminv(1-voxel_pval))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(fakecorrsz);
    % zero accounts for empty maps
    max_clust_info(permi) = max([ 0 ...
        cellfun(@numel,clustinfo.PixelIdxList) ]); 
    % using cellfun here eliminates the need for a slower loop over cells
end

% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
% find islands and remove those smaller than cluster size threshold
clustinfo = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
clust_threshold = prctile(max_clust_info,100-cluster_pval*100);

% identify clusters to remove
whichclusters2remove = find(clust_info<clust_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end


% plots with multiple comparisons corrections

% now find clusters in the real thresholded zmap
% if they are "too small" set them to zero
islands = bwconncomp(zmap);
for i=1:islands.NumObjects
    % if real clusters are too small, remove them by setting to zero!
    if numel(islands.PixelIdxList{i}==i)<clust_threshold
        zmap(islands.PixelIdxList{i})=0;
    end
end




figure
climdb  = [0 200];
contourf(times2save,frex,realmean, 20,'linecolor','none')
set(gca,'clim',climdb,'yscale','log','ytick',logspace(log10(min_frex),...
    log10(max_frex),6),'yticklabel',round(logspace(log10(min_frex),...
    log10(max_frex),6)*10)/10,'xlim', [-500 2000])
axis square
title('Averaged ITPCz')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
hold on
contour(times2save,frex,zmapthresh, 1,'linecolor','k')
colorbar
ylabel(colorbar, 'ITPCz')
cbh.FontSize = 28
% ylabel(colorbar, 'dB  from baseline', 'FontSize', 18)
ax = gca;
ax.FontSize = 20;
x0=10;
y0=10;
width=800;
height=500;
print('figure_itpcz_reflected_rmtl', '-dpng', '-r1000')

%% Theta
% Definition of TF windows

time_windows = [ 300 580 ]; % 3 time windows
freq_windows = [ frex(17) frex(20) ]; % 3 frequency windows

% find indices corresponding to time and frequency windows

timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));

for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(times2save-time_windows(i,j)));
    end
end

%for the frequency windows
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(frex-freq_windows(i,j)));
    end
end

for condi=1:size(tf,2) %
    
    % pointer to stats file
    statsfilename=(['statistics_file_TF_itpc_theta_reflected_rmtl_', ...
        num2str(condi),'.txt']);
    
    fid=fopen(statsfilename,'w');
    
    for subno=0:size(tf,1)
        
        
        % Write out column variable name or subject number.
        if subno==0
            fprintf(fid,'subnum\t');
        else
            fprintf(fid,'%g\t',subno);
        end
        
        for ti=1:size(time_windows,1)
            
            for fi=1:size(freq_windows,1)
                
                % Write out the column variable name of data.
                if subno==0
                    
                    fprintf(fid,[ num2str(time_windows(ti,1)) ...
                        num2str(time_windows(ti,2)) '_f ' ...
                        num2str(freq_windows(fi,1)) ...
                        num2str(freq_windows(fi,2)) '\t' ]);
                    
                else
                    
                    % get data from large window
                    % Note: This is condition-specific; you could make it
                    % condition-average by taking the mean across the 
                    % 2nd dimension.
                    tempdat = squeeze(itpcz(subno,condi,1,...
                        freqidx(fi,1):freqidx(fi,2),...
                        timeidx(ti,1):timeidx(ti,2)));
                    
                    % find max point (must first vectorize otherwise 
                    % you get max for each line)
                    [junk,maxpoint]=max(tempdat(:));
                    
                    % convert that back to 2d coordinates
                    [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                    
                    % get the indices from the larger matrix, not the 
                    % indices of the smaller TF window.
                    maxF = maxF+freqidx(1)-1;
                    maxT = maxT+timeidx(1)-1;
                    
                    % Now write out data and max time/frequency to file.
                    % 6 so as to have a proportional number of frequency 
                    % (3) and timepoints in the subject TF window
                    fprintf(fid,'%g\t',mean(mean(tempdat))); 
                end
            end
        end
        
        % CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
        fprintf(fid,'\n'); 
        if subno > 0
            maxfreqtheta(subno)=frex(maxF);
        else
            continue
        end
    end % end subject loop
    
    fclose(fid);
end % end condition loop

maxfreqtheta=maxfreqtheta';

save('maxFreq_theta.txt', 'maxfreqtheta', '-ascii', '-double', '-tabs')


%% Delta
% Definition of TF windows

time_windows = [ 370 595 ]; % 3 time windows
freq_windows = [ frex(10) frex(14) ]; % 3 frequency windows

% find indices corresponding to time and frequency windows

timeidx  = zeros(size(time_windows));
freqidx  = zeros(size(freq_windows));

for i=1:size(time_windows,1)
    for j=1:2
        [~,timeidx(i,j)] = min(abs(times2save-time_windows(i,j)));
    end
end

%for the frequency windows
for i=1:size(freq_windows,1)
    for j=1:2
        [~,freqidx(i,j)] = min(abs(frex-freq_windows(i,j)));
    end
end

for condi=1:size(tf,2) %
    
    % pointer to stats file
    statsfilename=(['statistics_file_TF_itpc_delta_reflected_rmtl',...
        num2str(condi), '.txt']);
    
    fid=fopen(statsfilename,'w');
    
    for subno=0:size(tf,1)
        
        
        % Write out column variable name or subject number.
        if subno==0
            fprintf(fid,'subnum\t');
        else
            fprintf(fid,'%g\t',subno);
        end
        
        for ti=1:size(time_windows,1)
            
            for fi=1:size(freq_windows,1)
                
                % Write out the column variable name of data.
                if subno==0
                    
                    fprintf(fid,[ num2str(time_windows(ti,1)) ...
                        num2str(time_windows(ti,2)) '_f ' ...
                        num2str(freq_windows(fi,1)) ...
                        num2str(freq_windows(fi,2)) '\t' ]);
                    
                else
                    
                    % get data from large window
                    % Note: This is condition-specific; 
                    % you could make it condition-average by taking the 
                    % mean across the 2nd dimension.
                    tempdat = squeeze(itpcz(subno,condi,1,freqidx(fi,1):...
                        freqidx(fi,2),timeidx(ti,1):timeidx(ti,2)));
                    
                    % find max point (must first vectorize otherwise 
                    % you get max for each line)
                    [junk,maxpoint]=max(tempdat(:));
                    
                    % convert that back to 2d coordinates
                    [maxF,maxT]=ind2sub(size(tempdat),maxpoint);
                    
                    % get the indices from the larger matrix, not the 
                    % indices of the smaller TF window.
                    maxF = maxF+freqidx(1)-1;
                    maxT = maxT+timeidx(1)-1;
                    
                    % Now write out data and max time/frequency to file.
                    % 6 so as to have a proportional number of frequency 
                    %(3) and timepoints in the subject TF window
                    fprintf(fid,'%g\t',mean(mean(tempdat))); 
                end
            end
        end
        
        % CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
        fprintf(fid,'\n');
        if subno > 0
            maxfreqdelta(subno)=frex(maxF);
        else
            continue
        end
    end % end subject loop
    
    fclose(fid);
end % end condition loop

maxfreqdelta=maxfreqdelta';

save('maxFreq_delta.txt', 'maxfreqdelta', '-ascii', '-double', '-tabs')

%% End