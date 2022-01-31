%% Load data

matfiles = dir(fullfile(['*.mat']));
numfiles = length(matfiles);

for subno = 1:numfiles
    
    load(matfiles(subno).name)
    
    % initialize output time-frequency data
    tf = zeros(3,num_frex,size(times2save,2),size(LFPemo.data, 3)); % 6 for the number of conditions
    
    
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
    s = logspace(log10(nCycRange(1)),log10(nCycRange(2)),num_frex)./(2*pi.*frex);
    t = -2:1/LFPemo.srate:2;
    halfwave = floor((length(t)-1)/2);
    
    % Data convolution parameters
    nData = size(LFPemo.data, 3)*LFPemo.pnts;
    nWave = length(t);
    nConv = nData+nWave-1;
    
    % Data wavelets
    cmwX = zeros(num_frex,nConv);
    for fi=1:num_frex
        cmw = fft(  exp(1i*2*pi*frex(fi).*t) .* exp( (-t.^2)/(2*s(fi)^2) )  ,nConv);
        cmwX(fi,:) = cmw ./ max(cmw); % amplitude-normalize in the frequency domain
    end
    
    % indices for the baseline
    baseidx = dsearchn(LFPemo.times',basetime');
    
    % Careful, do this in a condition-specific way
    for chani = 1:size(LFPemo.data,1)
        dataX = fft( reshape(LFPemo.data(chani,:,:),1,nData) ,nConv); %% concatenates all the trials of the first channel
        
        for fi=1:num_frex
            as = ifft( dataX.*cmwX(fi,:) );
            as = as(halfwave+1:end-halfwave);
            as = reshape(as,1,LFPemo.pnts, size(LFPemo.data, 3));
            
            % condition-average baseline
            basepow = mean(squeeze(mean(abs(as(:,baseidx(1):baseidx(2),:)).^2,2)),1);
            %baseitpc = mean(abs(mean(squeeze(exp(1i*angle(as_base))),2)),1);
            
            
            for triali = 1:size(LFPemo.data,3)
                tf(chani,fi,:,triali) = 10*log10( (abs(as(:,times2saveidx,triali)).^2) ./ basepow ); %
            end % end trial loop
        end % end frequencies loop
    end% end channel loop
    tfstruct(subno).data = tf;
    disp(['subject ', num2str(subno), ' is done'])
end% subject loop



%% Extract power from the tf windows of interest

%% Power

%% Delta 1
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
for subno=1:size(tfstruct,2)
    
    load(matfiles(subno).name)
    
    tf = tfstruct(subno).data;
    tempdat = zeros(size(tf,4), 3);
    for triali=1:size(tf,4)
        tempdat(triali, 1) = subno;
        tempdat(triali, 2) = mean(mean(tf(1,freqidx(1):freqidx(2),timeidx(1):timeidx(2), triali)));
        tempdat(triali, 3) = LFPemo.condition(triali);
        
    end
    tfstruct(subno).meantfdat = tempdat;
end


    % pointer to stats file
statsfilename=(['statistics_file_TF_delta_reflected_rmtl.txt']);

fid=fopen(statsfilename,'w');

for subno=1:size(tfstruct,2)
    
        load(matfiles(subno).name);
        tempdat = tfstruct(subno).meantfdat;
        for triali=1:size(tempdat, 1)
            fprintf(fid,'%g\t', tempdat(triali,:)); % 6 so as to have a proportional number of frequency (3) and timepoints in the subject TF window
            fprintf(fid,'\n'); %% CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
            
        end
    end
    
    

    fclose(fid);


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
for subno=1:size(tfstruct,2)
    
    load(matfiles(subno).name)
    
    tf = tfstruct(subno).data;
    tempdat = zeros(size(tf,4), 3);
    for triali=1:size(tf,4)
        tempdat(triali, 1) = subno;
        tempdat(triali, 2) = mean(mean(tf(1,freqidx(1):freqidx(2),timeidx(1):timeidx(2), triali)));
        tempdat(triali, 3) = LFPemo.condition(triali);
        
    end
    tfstruct(subno).meantfdat = tempdat;
    tfstruct(subno).meantfdat_delta2 = tempdat;

end


    % pointer to stats file
statsfilename=(['statistics_file_TF_alpha_reflected_rmtl.txt']);

fid=fopen(statsfilename,'w');

for subno=1:size(tfstruct,2)
    
        load(matfiles(subno).name);
        tempdat = tfstruct(subno).meantfdat;
        for triali=1:size(tempdat, 1)
            fprintf(fid,'%g\t', tempdat(triali,:)); % 6 so as to have a proportional number of frequency (3) and timepoints in the subject TF window
            fprintf(fid,'\n'); %% CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
            
        end
    end
    
    

    fclose(fid);


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
for subno=1:size(tfstruct,2)
    
    load(matfiles(subno).name)
    
    tf = tfstruct(subno).data;
    tempdat = zeros(size(tf,4), 3);
    for triali=1:size(tf,4)
        tempdat(triali, 1) = subno;
        tempdat(triali, 2) = mean(mean(tf(1,freqidx(1):freqidx(2),timeidx(1):timeidx(2), triali)));
        tempdat(triali, 3) = LFPemo.condition(triali);
        
    end
    tfstruct(subno).meantfdat = tempdat;
    tfstruct(subno).meantfdat_delta3 = tempdat;

end


    % pointer to stats file
statsfilename=(['statistics_file_TF_beta_reflected_rmtl.txt']);

fid=fopen(statsfilename,'w');

for subno=1:size(tfstruct,2)
    
        load(matfiles(subno).name);
        tempdat = tfstruct(subno).meantfdat;
        for triali=1:size(tempdat, 1)
            fprintf(fid,'%g\t', tempdat(triali,:)); % 6 so as to have a proportional number of frequency (3) and timepoints in the subject TF window
            fprintf(fid,'\n'); %% CAREFUL THIS IS CRITICAL FOR GOOD ORGANIZATIONS OF LINES !!!
            
        end
    end
    
    

    fclose(fid);

