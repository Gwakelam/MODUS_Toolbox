%% Process_data_Test_Set
%
% Processes multiple datafile outputs from CUBE query script in target folder,
% can specify target variables to include & process.
% Also finds Hm0 data from CUBElite and adds to dataset
%
%
% 19/04/18      First Version
% 15/05/18      Revised Version (changes not recorded...) GW
% 24/05/18      Adapted for use with new processing tools & multiple TKE
%               density methods
% 07/06/18      GW - v2 - Added handling of changing cell widths & blanking distances 
%               removed plotting options for ADCP data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear

tide        = 'flood';
direction   = {'X'};    % Only process 1 direction at a time for SBD
orientation = {'FF'};   % Turbine orientation - can process multiple orientations in order - must match name of sub-folder(s)
vars        = {'B1'};   % variable name for measured data must be specified

% Data input/load
% Flood
if strcmp(tide, 'flood')
    FolderIn = 'C:\3. ReDAPT Matlab Work\Higher Order No Wave Paper\Lscale - xyz test\Query_xyz_no_vbins\flood_x_all_vel';
    FolderOut = 'C:\3. ReDAPT Matlab Work\Higher Order No Wave Paper\Lscale - xyz test\Query_xyz_no_vbins\flood_x_processed_v2';
    % Ebb
elseif strcmp(tide, 'ebb')
    FolderIn = '';
    FolderOut = '';
    % Test
elseif strcmp(tide, 'test_flood')
    FolderIn = 'C:\3. ReDAPT Matlab Work\Higher Order No Wave Paper\Lscale - xyz test\Query_x_06_06_18\flood_2_1_ms\FF\vels001\1';
    FolderOut = 'C:\3. ReDAPT Matlab Work\Higher Order No Wave Paper\Lscale - xyz test\Query_x_06_06_18\Processed';
elseif strcmp(tide, 'test_ebb')
    FolderIn = '';
    FolderOut = '';
end

% Path to processing scripts
scriptloc = 'C:\3. ReDAPT Matlab Work\Own_Codes';
addpath(genpath(scriptloc));

% Switches
switch_result_only = 0;
switch_plot = 1;

% Define directions/beams to be processed
% Function switch for each direction/beam
LSswitch    = 1;   % Lengthscale
TIswitch    = 1;   % TI
Dissswitch  = 0;   % TKE production (structure function)
LStswitch   = 0;   % Transform of ADCP along-beam lengthscales to flow axes
AmpSwitch   = 1;   % Extract beam amplitudes
CorrSwitch  = 1;   % Extract beam correlation
QC_on       = 0;   % QC interp switch

% Set upper limit for slack water
u_slack = 0.5;  %m/s
options = [];   % Empty options struct - uses defult vals for all funcs


for o = 1:length(orientation)  % Orientation loop
    
    direct = strcat(FolderIn, '\', orientation{o});
    cd(direct);
    files = dir;
    files = files(3:end);
    nfile = length(files);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% MAIN LOOP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1:nfile  % Loop through each data file
        
        fprintf('[INFO] - Loading datafile %d of %d \n', j, nfile);
        load(files(j).name)
        sensors = fields(Data.(direction{1}));
        
        for s = 1:length(sensors)
            fprintf(['[INFO] - Processing sensor ' sensors{s} '\n']);
  
            % Check sample rates
            Fs_chk = unique(Data.(direction{1}).(sensors{s}).CF(:,6));
            Fs_chk =  Fs_chk(~isnan( Fs_chk));   % Ignore occasional NaN's in sample rate field
            
            for f = 1:length(Fs_chk)   % Loop through each sample rate
                fprintf('[INFO] - Processing %d Hz data \n', Fs_chk(f))
                
                clear Results Data_out
                
                % set index for same sample rate - not always conigious
                Fsind = find(Data.(direction{1}).(sensors{s}).CF(:,6) == Fs_chk(f));
                
%                 if length(Data_out.ref.Fs) ~= 1
%                     error('more than one sample rate in data')
%                 end
                
                %% GET REFVALS
                
                % test = Data.(direction{1}).(sensors{s}).Timestamp(Fsind);
                
                % Ref values
                Data_out.ref.Fs         = Fs_chk(f);
                Data_out.ref.orientation = orientation;
                Data_out.ref.Tref   = Data.Info.Turbine_Status.Timestamp;
                Data_out.ref.Uref   = Data.Info.Turbine_Status.ReferenceVelocity;
                Data_out.ref.Tstat  = (Data_out.ref.Tref(2) - Data_out.ref.Tref(1))*(24*60^2);  % Tstat length in sec
                Data_out.ref.timestamp = Tstat_reshape_v1(Data.(direction{1}).(sensors{s}).Timestamp(Fsind), ...
                                        Data_out.ref.Fs, Data_out.ref.Tstat);
                

                % Trim Uref & Tref to match actual dataset
                [Data_out.ref.Tref, Data_out.ref.Uref] = ...
                    trim_refdata_v2(Data_out.ref.Tref, Data_out.ref.Uref, Data_out.ref.timestamp, round(Data_out.ref.Tstat/60));
                
                %test
                switch_trimtest = 0;
                if switch_trimtest == 1
                    Tref = Data_out.ref.Tref;
                    Uref = Data_out.ref.Uref;
                    timestamp = Data_out.ref.timestamp;
                    hm0 = Data_out.ref.hm0;
                    Tstat = round(Data_out.ref.Tstat/60);
                end
                
                % Get hm0 values to match dataset
                Data_out.ref.hm0 = get_hm0_v1(Data_out.ref.Tref);
                
                if switch_plot == 1
                    figure
                    plot(Data.(direction{1}).(sensors{s}).Timestamp, Data.(direction{1}).(sensors{s}).B1(:, 20), Data_out.ref.Tref, Data_out.ref.Uref, ...
                        Data_out.ref.Tref, Data_out.ref.hm0)
                end
          
                % Get distance/depth bins - preserve all values as they can change during deployment
                Data_out.ref.first_cell = unique(Data.(direction{1}).(sensors{s}).CF(Fsind,5));
                Data_out.ref.first_cell_all = Tstat_reshape_v1(Data.(direction{1}).(sensors{s}).CF(Fsind,5), ...
                    Data_out.ref.Fs, Data_out.ref.Tstat);
        
                Data_out.ref.cell_size = unique(Data.(direction{1}).(sensors{s}).CF(Fsind,4));
                Data_out.ref.cell_size_all = Tstat_reshape_v1(Data.(direction{1}).(sensors{s}).CF(Fsind,4), ...
                    Data_out.ref.Fs, Data_out.ref.Tstat);

                nz = size(Data.(direction{1}).(sensors{s}).B1, 2);
                
                % Generate ensemble of bin distances, kill any ensemble where bin properties change during Tstat
                [Data_out.ref.zbin, Data_out.ref.binQC] = QC_cellbins_v1(Data_out.ref.cell_size_all, Data_out.ref.first_cell_all, nz);  
                
                if AmpSwitch
                    Data_out.Amp1.measured = Tstat_reshape_v1(Data.(direction{1}).(sensors{s}).Amp1, ...
                        Data_out.ref.Fs, Data_out.ref.Tstat);
                    Data_out.Amp1.mean = squeeze(nanmean(Data_out.Amp1.measured, 1));
                end
                
                if CorrSwitch
                    Data_out.Cor1.measured = Tstat_reshape_v1(Data.(direction{1}).(sensors{s}).Cor1, ...
                        Data_out.ref.Fs, Data_out.ref.Tstat);
                    Data_out.Cor1.mean = squeeze(nanmean(Data_out.Cor1.measured, 1));
                end
                
                
                %% SET VELOCITY BINS, INDEX TSTAT ENSEMBLES, Hm0 limit
                Data_out.ref.v_edges = [0.3 0.7 1.1 1.5 1.9 2.3 2.7 3.1 3.5 3.9];
                Data_out.ref.v_mids  = [0.5 0.9 1.3 1.7 2.1 2.5 2.9 3.3 3.7];
                
                Data_out.ref.vbins = discretize(abs(Data_out.ref.Uref), Data_out.ref.v_edges);
                
                
                %% PROCESS DATA
                
                for i = 1:length(vars)  % Loop through specified direction/beam variables
                    
                    fprintf('[INFO] - Processing %s data.\n', vars{i})
                    
                    % Reshape dataset
                    Data_out.(vars{i}).measured   = Tstat_reshape_v1(Data.(direction{1}).(sensors{s}).(vars{i}), ...
                        Data_out.ref.Fs, Data_out.ref.Tstat);
                    
                    % NaN treatment
                    badlim = 5; % max 5% Interpolated values
                    maxgap = 1; % Max 1 sample gap remaining after interpolation
                    
                    if QC_on
                        [Data_out.(vars{i}).measured, Data_out.(vars{i}).fail_QC]  = ...
                            interp_QC_v1(Data_out.(vars{i}).measured, badlim, maxgap);
                    end
                    
                    % Detrending (linear)
                    [Data_out.(vars{i}).detrended, Data_out.(vars{i}).mean, Data_out.(vars{i}).stdev] = detrend_lin3D_v2(...
                        Data_out.(vars{i}).measured, 1);
                    
                    [ns ne nz] = size(Data_out.(vars{i}).detrended);
                    
                    
                    % Pass only slack water Tstats to noise correction function (see prev work)
                    slacks = Data_out.(vars{i}).detrended(:, abs(Data_out.ref.Uref) <= u_slack, :);
                    
                    [Data_out.(vars{i}).noise_corr, Data_out.(vars{i}).f, Data_out.(vars{i}).PSD_slacks] = PSD_noisebias_v2(...
                        slacks, Data_out.ref.Fs, []);
                    
                    if switch_plot == 2
                        figure
                        loglog(Data_out.(vars{i}).f, Data_out.(vars{i}).PSD_slacks(:,1,20))
                    end
                    
                    % Calc integral lengthscales if enabled by LSswitch
                    if LSswitch(i)
                        Data_out.(vars{i}).Lscale_ACF = calc_Lscale_v1(Data_out.(vars{i}).detrended, Data_out.(vars{i}).mean, ...
                            Data_out.ref.Fs, ceil(0.9*ns), 0);
                        
                        [Data_out.(vars{i}).Lscale_Macro, ] = calc_macro_Lscale_v1(Data_out.(vars{i}).detrended, ...
                            Data_out.(vars{i}).mean,  Data_out.ref.Fs, 0.02);
                    end
                    
                    % Calc TI - Note that noise_corr output of PSD_noisebias_ is VARIANCE - must be square rooted to give
                    % std dev input for function. Single value correction used here -
                    % no difference for different depths and Tstat ensembles
                    corr_range = [10:30];
                    if TIswitch(i)
                        Data_out.(vars{i}).TI     = calc_TI_v2(Data_out.(vars{i}).stdev, ...
                            sqrt(nanmean(Data_out.(vars{i}).noise_corr(corr_range))), Data_out.(vars{i}).mean, []) ;
                    end
                    
                    if Dissswitch(i)
                        Data_out.(vars{i}).TKE_Diss = calc_Diss_v1(Data_out.(vars{i}).detrended, Data_out.ref.zbin, options);
                    end
                    
                    
                end
                
       
                
                %% Save file
                if ~isdir(FolderOut)
                    mkdir(FolderOut)
                end
                
                namestr = strcat(sensors{s}, '_OUT_', direction{1}, '_', orientation{o}, '_', num2str(Data_out.ref.Fs), 'Hz' );
                cd(FolderOut)
                save(namestr, 'Data_out')
                
            end
            
        end
        
        %% Test Plots
        %z = Data_out.ref.zbin;
        
        if switch_plot == 9
 
            
            %% Plot Beam velocities
            figure,clf
            pcolor(Data_out.B1.mean), shading flat
            colorbar
            caxis(1*[-3 3])
            ylabel('Beam 1')
            title('veloctiy [m/s]')
           

            
            %% Plot Amplitudes
            lim = 50;
            figure,clf
            pcolor(Data_out.Amp1.mean), shading flat
            colorbar
            caxis(lim*[-0 1])
            ylabel('Beam 1')
            title('Amplitude')
            
            %% Plot correlation
            lim = 90;
            figure,clf
            pcolor(Data_out.Cor1.mean), shading flat
            colorbar
            caxis(lim*[-0 1])
            ylabel('Beam 1')
            title('Correlation')
            
            %% Plot TI
            figure,clf
            pcolor(Data_out.B1.TI), shading flat
            colorbar
            caxis([-0 15])
            ylabel('Beam 1')
            title('TI')      
            
             %% Plot Lengthscale
            figure,clf
            pcolor(abs(Data_out.B1.Lscale_ACF)), shading flat
            colorbar
            caxis([-0 20])
            ylabel('Beam 1')
            title('LscaleACF')   
            
        end
        
        
        
    end
    
end

%% Test Plots
debug = 0;

if debug
    figure,clf
    pcolor(Data.(direction{1}).(sensors{s}).Amp1' ), shading flat
    colorbar
    caxis(120*[0 1])
    ylabel('Depth bin')
    title('Beam 2 Amplitude')
    
    figure,clf
    pcolor(Data.(direction{1}).(sensors{s}).B1'), shading flat
    colorbar
    caxis(1*[-1 1])
    ylabel('Depth bin')
    title('Beam2 Velocity')
    
    figure,clf
    pcolor(Data.(direction{1}).(sensors{s}).U'), shading flat
    colorbar
    caxis(3*[-1 1])
    ylabel('Depth bin')
    title('Transformed U Velocity')
    
    figure,clf
    pcolor((Data_out.B1.stdev.^2)'), shading flat
    colorbar
    caxis(0.02*[-1 1])
    ylabel('Depth bin')
    title('B2 Variance')
    
    figure,clf
    pcolor(Data_out.B1.mean'), shading flat
    colorbar
    caxis(0.8*[-1 1])
    ylabel('Depth bin')
    title('B2 Variance')
    

    
    figure,clf
    pcolor(Data_out.B1.Lscale_ACF'), shading flat
    colorbar
    caxis(4*[-1 1])
    ylabel('Depth bin')
    title('B2 Lengthscale') 
    
end


   