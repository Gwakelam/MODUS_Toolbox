% Takes processed CUBE SBD data (from Process_data_vX) and consolidates subsets
% of different results files based of flow velocity and hm0
%
%
%  OUTPUTS:
%
% 27/04/18      First Version
% 24/05/18      Adapted for use with test data set - GW

%% SETUP
clear
close all


tide = 'test_flood';

% Data input/load
% Flood 0
if strcmp(tide, 'flood_')
    FolderIn = '';
    FolderOut = '';
    % Ebb
elseif strcmp(tide, 'ebb')
    FolderIn = '';
    FolderOut = '';
    % Test
elseif strcmp(tide, 'test_flood')
    FolderIn = 'C:\3. ReDAPT Matlab Work\Higher Order No Wave Paper\Lscale - xyz test\Query_xyz_no_vbins\flood_x_processed_v2';
    FolderOut = 'C:\3. ReDAPT Matlab Work\Higher Order No Wave Paper\Lscale - xyz test\Query_xyz_no_vbins\flood_x_results_v2';
elseif strcmp(tide, 'test_ebb')
    FolderIn = ' xxx ';
    FolderOut = 'xxx ';
end

% Path to processing scripts
scriptloc = 'C:\3. ReDAPT Matlab Work\Own_Codes';
addpath(genpath(scriptloc));
% Supporting functions
funcloc = 'C:\3. ReDAPT Matlab Work\useful_funcs';
addpath(genpath(funcloc));

% Switches
switch_result_only = 0;
switch_plot = 0;

%% Set Subsets

cd(FolderIn);
files = dir('*OUT*');

 load(files(1).name)
 vmids = Data_out.ref.v_mids;
 vedge = Data_out.ref.v_edges;
 
 
 for i = 1:length(vmids)
     vname{i} = strcat('vbin_', num2str(i));
 end

hm0_max = 1.0;
hm0_min = 0;


%% Get Datafiles

  % Initialise empty results arrays
   for j = 1:length(vmids)
    
       Results.(vname{j}).ref.Fs     = [];
       
        Results.(vname{j}).ref.Tref  = [];
        Results.(vname{j}).ref.Uref  = [];
        Results.(vname{j}).ref.hm0   = [];
        
        Results.(vname{j}).ref.zbin     = [];  
        Results.(vname{j}).ref.zbinQC   = [];
        
        Results.(vname{j}).B1.mean      = [];
        Results.(vname{j}).B1.std       = [];
        Results.(vname{j}).B1.TI            = [];
        Results.(vname{j}).B1.noise_corr    = [];
        Results.(vname{j}).B1.f             = [];
        Results.(vname{j}).B1.Slack_PSD     = [];
        Results.(vname{j}).B1.Amp1mean      = [];
        Results.(vname{j}).B1.Cor1mean      = [];
        
        Results.(vname{j}).LS.B1_ACF    = [];
        
        Results.(vname{j}).LS.B1_Macro  = [];
        
 
    
    end

for i = 1:length(files)
    
    fprintf('[INFO] - Loading file %d out of %d \n', i, length(files))
    load(files(i).name)
    fprintf('[INFO] - Loaded, concatenating data\n')
    % Hm0 already dealt with in Query for now
    % hm0 = find((Data_out.ref.hm0 >= hm0_min) & (Data_out.ref.hm0 <= hm0_max));
    
    % Concat data
    for j = 1:length(vmids)
        
        vbins = Data_out.ref.vbins;
        vel = find(vbins == j);
        
        ind = vel; %hm0(ismember(hm0, vel));
        indtot(i,j) = length(ind); 
        
        %Results.(vname{j}).ref.Fs    = horzcat(Results.(vname{j}).ref.Fs, repmat(Data_out.ref.Fs, 1, indtot));
        %Results.(vname{j}).ref.orientation = horzcat(Results.(vname{j}).ref.orientation, repmat(Data_out.ref.orientation, 1, indtot));
        
        Results.(vname{j}).ref.Tref  = vertcat(Results.(vname{j}).ref.Tref, Data_out.ref.Tref(ind));
        Results.(vname{j}).ref.Uref  = vertcat(Results.(vname{j}).ref.Uref, Data_out.ref.Uref(ind));
        Results.(vname{j}).ref.hm0  = vertcat(Results.(vname{j}).ref.hm0, Data_out.ref.hm0(ind));
        
        Results.(vname{j}).B1.mean = padconcatenation(Results.(vname{j}).B1.mean, Data_out.B1.mean(ind, :), 1);
        Results.(vname{j}).B1.std = padconcatenation(Results.(vname{j}).B1.std, Data_out.B1.stdev(ind, :), 1);
        Results.(vname{j}).B1.TI = padconcatenation(Results.(vname{j}).B1.TI, Data_out.B1.TI(ind,:), 1);
  
        % Add PSD spectra / mean>?
        Results.(vname{j}).B1.f = padconcatenation(Results.(vname{j}).B1.f, Data_out.B1.f, 2);
        slack_PSD = nanmean(nanmean(Data_out.B1.PSD_slacks(:,:,10:20), 2), 3);
        Results.(vname{j}).B1.Slack_PSD = padconcatenation(Results.(vname{j}).B1.Slack_PSD, slack_PSD, 2);
        Results.(vname{j}).B1.noise_corr = padconcatenation( Results.(vname{j}).B1.noise_corr, Data_out.B1.noise_corr, 1);
        
        % Amp & Corr
        Results.(vname{j}).B1.Amp1mean = padconcatenation(Results.(vname{j}).B1.Amp1mean, Data_out.Amp1.mean(ind, :), 1);
        Results.(vname{j}).B1.Cor1mean = padconcatenation(Results.(vname{j}).B1.Cor1mean, Data_out.Cor1.mean(ind, :), 1);
        
        % Lengthscales
        Results.(vname{j}).LS.B1_ACF = padconcatenation(Results.(vname{j}).LS.B1_ACF, Data_out.B1.Lscale_ACF(ind, :), 1);
        Results.(vname{j}).LS.B1_Macro = padconcatenation(Results.(vname{j}).LS.B1_Macro, Data_out.B1.Lscale_Macro(ind, :), 1);
       

    end
    
    Results.ref_orig = Data_out.ref;
    Results.counts = indtot;
    
end

if strcmp(tide, 'flood_01')
    Results_Flood1 = Results;
    if ~isdir(FolderOut)
        mkdir(FolderOut)
    end
    cd(FolderOut)
    save('Results_Flood1', 'Results_Flood1')
    
    
elseif strcmp(tide, 'ebb')
    Results_Ebb = Results;
    if ~isdir(FolderOut)
        mkdir(FolderOut)
    end
    cd(FolderOut)
    save('Results_Ebb', 'Results_Ebb')
    
    elseif strcmp(tide, 'test_flood')
    Results_Flood_test = Results;
    if ~isdir(FolderOut)
        mkdir(FolderOut)
    end
    cd(FolderOut)
    save('Results_Flood_test', 'Results_Flood_test')
    
    elseif strcmp(tide, 'test_ebb')
    Results_Ebb = Results;
    if ~isdir(FolderOut)
        mkdir(FolderOut)
    end
    cd(FolderOut)
    save('Results_Ebb_test', 'Results_Ebb_test')
    
end
    






