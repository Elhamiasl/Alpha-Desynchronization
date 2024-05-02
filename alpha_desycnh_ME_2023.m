%% **Alpha Desynchronization in 6-, 9-, and 12-months-old Infants**
    % After pre-processing the data and obtaining FFT and desynchronization values per electrode and frequency for each subject,
        %Use this code
         % to generate spectra plot to visualize the pre-processed EEG data
            % for cross-sectional and longitudinal sample
         % to generate  headplots illustrating the distribution of alpha
             % desynchronization across whole brain
         % to run one sample t-test against the null hypothesis value of zero
             %(0) to examine significant desynchronization across each electrode per group and per
             %condition                   
%================================================ ME, November-December 2023
%% Amplitude Spectrum plot

% clear workspace and pull epoched FFT data from pre-processed data and average over trials of each condition per age group
    clear
    
    for grp = 1:3 % define group names based on the age specified in file names
        if grp==1, gID = '6'; elseif grp==2, gID = '9'; elseif grp==3, gID = '12'; end
        data = getfilesindir(pwd, ['*_' gID '_*.mat']); % get .mat files

        face_data_tmp = []; % Create empty output arrays for each condition  to prevent errors and avoid overwriting existing outputs
        obj_data_tmp = [];
     
        for sub = 1:size(data,1) 
            matName = deblank(data(sub,:)); 
                load(matName, 'reseg_1sFFT') % % load the current subject's 1-sec epoched FFT datamat
                if isfield(reseg_1sFFT,'fa06') == 1 % if the current subject has data for this condition (here Face (fa06))
                    face_data_tmp(:,:,sub) = mean(reseg_1sFFT.fa06,3); end % average over trial/segments => elec x freq x subject 
                if isfield(reseg_1sFFT,'un05') == 1 
                    obj_data_tmp(:,:,sub) = mean(reseg_1sFFT.un05,3); end                 
        end
        gm_faces(:,:,grp) = mean(face_data_tmp,3); % average over subjects and get grand mean => elec x freq x age group 
        gm_objs(:,:,grp) = mean(obj_data_tmp,3);
    end
  

% plot amplitude spectra by frequency   
    faxis = 0:1/1:250 - 1/1; % Calculate frequency axis
    figure
    gm_face_6mo = gm_faces(:,:,1); % pull the required data
    plot(faxis(2:21),mean(mean(gm_face_6mo([63 64  66 67 68 73 74], 2:21,:),3),1),'b','linewidth' , 6) %amplitude spectra by 1-20 Hz averaged over occipital electrodes
    hold on
    gm_face_9mo = gm_faces(:,:,2);
    plot(faxis(2:21),mean(mean(gm_face_9mo([63 64  66 67 68 73 74], 2:21,:),3),1),'g','linewidth' , 6) 
    gm_face_12mo = gm_faces(:,:,3);
    plot(faxis(2:21),mean(mean(gm_face_12mo([63 64  66 67 68 73 74], 2:21,:),3),1),'r','linewidth' , 6) 
    %xline(5)
    %set(gca,'FontSize',26)
    title('Face')
    legend('6mo', '9mo', '12mo')
    
%% Alpha Desynchronization

% clear workspace and pull desynchronization data from pre-processed data and average over trials of each condition per age group
    clear

    face = cell(1,3) % creat cell array for each condition
    object = cell(1,3)

    for grp = 1:3 % define group names based on the age specified in file names
        if grp==1, gID = '6'; elseif grp==2, gID = '9'; elseif grp==3, gID = '12'; end
        data = getfilesindir(pwd, ['*_' gID '_*.mat']); % all data 

        face_desynch_tmp = []; %Create empty output arrays for each condition
        obj_desynch_tmp = [];
     
        for sub = 1:size(data,1)
            matName = deblank(data(sub,:));
        
                load(matName, 'desynchSpecData')
                if isfield(desynchSpecData,'fa06') == 1 % if the current subject has data for this condition (here Face condition(fa06))
                    face_desynch_tmp(:,:,sub) = mean(desynchSpecData.fa06,3); end % average over trial/segments => elec x freq x subject 
                if isfield(desynchSpecData,'un05') == 1
                    obj_desynch_tmp(:,:,sub) = mean(desynchSpecData.un05,3); end                 
        end
    
        % average data over desired frequency band (Here 6 to 9 Hz which corresponds to 7 to 10 frequency bins) 
        % make sure to compute faxis in order to select the bins correspond to specific frequencies
        face{grp} = squeeze(mean(face_desynch_tmp(:,7:10,:),2)); %  => elec  x participant in each age group 
        obj{grp} = squeeze(mean(obj_desynch_tmp(:,7:10,:),2));
        
        % average over subjects and desired frequency and get grand means for topoplots
        gm_desynch_faces_topo (:,:,grp) =mean(mean(face_desynch_tmp(:,7:10,:),2),3); %  => elec  x age group 
        gm_desynch_obj_topo (:,:,grp) =mean(mean(obj_desynch_tmp(:,7:10,:),2),3); 
    end
  
%% Desynchronization grand mean topographies

    load locsEEGLAB109HCL % # load 109 electrode montage
    topoplot(gm_desynch_faces_topo (:,:,1), locsEEGLAB109HCL)  % modify the maps lo/hi limits if required
    title('desynch-FACE-6mo')
    colorbar

    topoplot(gm_desynch_obj_topo (:,:,1), locsEEGLAB109HCL)  
    title('desynch-OBJECT-6mo')
    colorbar
%% One sample t-test against the null hypothesis value of zero (0) across each electrode to examine significant desycnhronizations.

% for each condition and age group, loop over each electrode and run t-test vs 0 
    for elec = 1:109
        [~,pvalue(elec,1),~,stats] = ttest(face{2}(elec, :)); 
        tvals(elec,1) = stats.tstat;
    end

% determine which electrodes have significant decrease in alpha power, at p < .05
    sig_chans_desynch = zeros(109,1); 
    for elec = 1:109
        if pvalue(elec) <=.05 && tvals(elec) < 0, sig_chans_desynch(elec) = 1; end 
    end
    
% plot significant pvalues
    load locsEEGLAB109HCL % load 109 electrode montage
    topoplot(sig_chans_desynch, locsEEGLAB109HCL, 'electrodes' ,'on', 'style', 'blank')
    title('Desynchronization')
    colorbar