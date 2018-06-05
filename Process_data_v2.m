%% Process_data_v2  
%
% Processes multiple datafile outputs from CUBE query script in target folder,
% can specify target variables to include & process.
% Also finds Hm0 data from CUBElite and adds to dataset
%
%
% 19/04/18      First Version
% 15/05/18      Revised Version (changes not recorded...) GW

%% SETUP
tide = 'ebb';

% Data input/load
% Flood
if strcmp(tide, 'flood')
    FolderIn = 'C:\3. ReDAPT Matlab Work\For Manuel\Data_03_Apr_Flood';
    FolderOut = 'C:\3. ReDAPT Matlab Work\For Manuel\Flood_Processed';
% Ebb
elseif strcmp(tide, 'ebb')
    FolderIn = 'C:\3. ReDAPT Matlab Work\For Manuel\Data_03_Apr_Ebb2';
    FolderOut = 'C:\3. ReDAPT Matlab Work\For Manuel\Ebb_Processed';
end

% Path to processing scripts
scriptloc = 'C:\3. ReDAPT Matlab Work\Own_Codes';
addpath(genpath(scriptloc));

% Switches
switch_result_only = 0;
switch_plot = 0;

% Define directions/beams to be processed
vars = {'B1', 'B2', 'B3', 'B4', 'U'};
% Switch for lengthscale/TI for each direction/beam
LSswitch = [1 1 1 1 1];     % Lengthscale
TIswitch = [0 0 0 0 1];     % TI
% Set upper limit for slack water
u_slack = 0.5;  %m/s

cd(FolderIn);
files = dir('*Data_*');
nfile = length(files);

for j = 1:nfile
    
    clear Results Data_out
    
    fprintf('[INFO] - Loading datafile %d of %d \n', j, nfile);
    cd(FolderIn);
    load(files(j).name)
    
    
    % Get Deployment name from data field
    dep = fields(Data.ADCP);
    
    % Check sample rates
    Data_out.ref.Fs = unique(Data.ADCP.(dep{1}).CF(:,13));
    
    if length(Data_out.ref.Fs) ~= 1
        error('more than one sample rate in data')
    end
    
    %% GET REFVALS
    
    % Ref values
    Data_out.ref.Tref = Data.Info.Turbine_Status.Timestamp;
    Data_out.ref.Uref = Data.Info.Turbine_Status.ReferenceVelocity;
    Data_out.ref.Tstat = (Data_out.ref.Tref(2) - Data_out.ref.Tref(1))*(24*60^2);  % Tstat length in sec
    Data_out.ref.timestamp = Tstat_reshape_v1(Data.ADCP.(dep{1}).Timestamp, ...
        Data_out.ref.Fs, Data_out.ref.Tstat);
    
    % Get hm0 values to match dataset
    Data_out.ref.hm0 = get_hm0_v1(Data_out.ref.Tref);
    
    % Trim Uref & Tref to match actual dataset
    [Data_out.ref.Tref, Data_out.ref.Uref, Data_out.ref.hm0 ] = trim_refdata_v1(Data_out.ref.Tref, Data_out.ref.Uref, Data_out.ref.timestamp, Data_out.ref.hm0 );
    
    figure
    plot(Data.ADCP.(dep{1}).Timestamp, Data.ADCP.(dep{1}).U(:, 20), Data_out.ref.Tref, Data_out.ref.Uref, ...
        Data_out.ref.Tref, Data_out.ref.hm0)
    
    % Install depth can vary over deployment, can be incorrect values that
    % are very low. This takes mean depth as a single value
    Data_out.ref.install_depth = nanmean(Data.ADCP.(dep{1}).CF( Data.ADCP.(dep{1}).CF(:,9) > 30 ,9));
    
    Data_out.ref.first_cell = unique(Data.ADCP.(dep{1}).CF(:,14)); % Does not account for frame height
    Data_out.ref.cell_size = unique(Data.ADCP.(dep{1}).CF(:,19));
    
    nz = size(Data.ADCP.(dep{1}).(vars{1}), 2);
    Data_out.ref.zbin = Data_out.ref.first_cell:Data_out.ref.cell_size: Data_out.ref.first_cell + (nz-1)*Data_out.ref.cell_size;
    
    %% SET VELOCITY BINS, INDEX TSTAT ENSEMBLES, Hm0 limit
    Data_out.ref.v_edges = [0 0.3 0.7 1.1 1.5 1.9 2.3 2.7 3.1];
    Data_out.ref.v_mids  = [0.15 0.5 0.9 1.3 1.7 2.1 2.5 2.9];
    
    Data_out.ref.vbins = discretize(abs(Data_out.ref.Uref), Data_out.ref.v_edges);
    
    
    %% PROCESS DATA
    
    for i = 1:length(vars)
        
        fprintf('[INFO] - Processing %s data.\n', vars{i})
        
        Data_out.(vars{i}).measured   = Tstat_reshape_v1(Data.ADCP.(dep{1}).(vars{i}), ...
            Data_out.ref.Fs, Data_out.ref.Tstat);
        
        
        [Data_out.(vars{i}).detrended, Data_out.(vars{i}).mean, Data_out.(vars{i}).stdev] = ...
             detrend_lin3D_v2(Data_out.(vars{i}).measured, 1);
        
        
        % Pass only slack water Tstats to noise correction function (see prev work)
        slacks = Data_out.(vars{i}).detrended(:, abs(Data_out.ref.Uref) <= u_slack, :);
        
        Data_out.(vars{i}).noise_corr = PSD_noisebias_v1(slacks, Data_out.ref.Fs, []);
        
        
        if LSswitch(i)
            Data_out.(vars{i}).Lscale = calc_Lscale_v1(Data_out.(vars{i}).detrended, Data_out.(vars{i}).mean, ...
                Data_out.ref.Fs, 125, 0);
        end
        
        % Calc TI - Note that noise_corr output of PSD_noisebias_ is VARIANCE - must be squared to give
        % std dev input for function
        if TIswitch(i)
            Data_out.(vars{i}).TI     = calc_TI_v2(Data_out.(vars{i}).stdev , nanmean(Data_out.(vars{i}).noise_corr(1:20)).^2, ...
                Data_out.(vars{i}).mean, []) ;
        end
        
    end
    
    
    %% Calc Rstress & TKE
    
    % Create input structs for
    
    std_vb.v1 =  Data_out.B1.stdev;
    std_vb.v2 =  Data_out.B2.stdev;
    std_vb.v3 =  Data_out.B3.stdev;
    std_vb.v4 =  Data_out.B4.stdev;
    
    % with depth varying noise correction - previous value interpolation used -
    % find better method!
    std_noise.v1 = interp1gap(Data_out.B1.noise_corr, 'previous');
    std_noise.v2 = interp1gap(Data_out.B2.noise_corr, 'previous');
    std_noise.v3 = interp1gap(Data_out.B3.noise_corr, 'previous');
    std_noise.v4 = interp1gap(Data_out.B4.noise_corr, 'previous');
    
    % use default options
    options = [];
    
    [Data_out.TKE, Data_out.Rstress] = calc_TKE_Rstress_v2(std_vb, std_noise, options);
    
    
    %% Save file
    if ~isdir(FolderOut)
        mkdir(FolderOut)
    end
    
    namestr = strcat(files(1).name(19:32),'_OUT_',  num2str(j));
    cd(FolderOut)
    save(namestr, 'Data_out') 

    
end
    %% Plotting
    
    if switch_plot
        
    v_bin = fields(Results);
    k = 2;
    
    figure
    plot(Results.(v_bin{k}).U, Results.ref.cellheight, ':.')
    legend('U mean vel',  'location', 'best')
    xlabel('mean inflow [m/s]'); ylabel('dist from bed [m]');
    xlim([0 3])
    grid on; set(gcf,'color','white');
    
    figure
    plot( Results.(v_bin{k}).TI, Results.ref.cellheight,':.')
    legend('TI_U',  'location', 'best')
    xlabel('TI [%]');  ylabel('dist from bed [m]');
    grid on; set(gcf,'color','white');
    xlim([0 30])
    
    figure
    plot(Results.(v_bin{k}).Rxz, Results.ref.cellheight, ':.', Results.(v_bin{k}).Ryz, Results.ref.cellheight, ':.')
    legend('Rxz', 'Ryz', 'location', 'best')
    xlabel('Reynolds Stress [(m/s)^2]');  ylabel('dist from bed [m]');
    grid on; set(gcf,'color','white');
    
    figure
    plot(Results.(v_bin{k}).Lscale_B1, Results.ref.cellheight, ':.', Results.(v_bin{k}).Lscale_B2, Results.ref.cellheight, ':.', ...
        Results.(v_bin{k}).Lscale_B3, Results.ref.cellheight, ':.', Results.(v_bin{k}).Lscale_B4, Results.ref.cellheight, ':.')
    legend('B1', 'B2', 'B3', 'B4',  'location', 'best')
    xlabel('along-beam integral lengthscale [m]');  ylabel('dist from bed [m]');
    grid on; set(gcf,'color','white');
    
    figure
    plot(Results.(v_bin{k}).TKE, Results.ref.cellheight, ':.')
    legend('TKE',  'location', 'best')
    xlabel('TKE'); ylabel('dist from bed [m]');
    xlim([0 0.35])
    grid on; set(gcf,'color','white');
    
    
    end