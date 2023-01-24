% working path: /home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Climate_variability_carbon_cycle/Analysis
%% ---- MATLAB NC Lib ---
clear global *;
clearvars;
sys_id='Linux_Dell_dlocal'; % Linux_Dell, ra_T420, tiberino

global data_src;
global lib_path;
global wrk_dir;

if strcmp(sys_id,'Linux_Dell')
    lib_path='/home/minchao/work/matlab_lib/';
    data_src='/media/minchao/DS1/DATA/obs/';
    wrk_dir='/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Climate_variability_carbon_cycle/Analysis/';
elseif strcmp(sys_id,'Linux_Dell_dlocal')
    lib_path='/home/minchao/work/matlab_lib/';
    data_src='/Data/obs';
    wrk_dir='/home/minchao/work/GoogleDriver/Research/a1_Studies_and_projects/Climate_variability_carbon_cycle/Analysis/';
elseif strcmp(sys_id,'ra_T420')
    lib_path='/home/minchao/NewWork/matlab_lib/';
    data_src='/media/minchao/DS1/DATA/obs/';
    wrk_dir='/home/minchao/Google_Drive/minchaowu.acd@gmail.com/Research/a1_Studies_and_projects/Climate_variability_carbon_cycle/Analysis';
elseif strcmp(sys_id,'tiberino')
    lib_path='/home/minchao/work/matlab_lib/';
    data_src='/data/fyris/minchao/';
    wrk_dir='/home/minchao/work/carbon_cycle_climate_variability/Analysis';
end

addpath(lib_path);
addpath(sprintf('%s/sm_grini/nctoolbox/nctoolbox-1.1.1/',lib_path));
setup_nctoolbox;
addpath(sprintf('%s/CDT-Chad_Greene/v1.01/cdt/',lib_path));
addpath(sprintf('%s/CDT-Chad_Greene/v1.01/cdt/cdt_data/',lib_path));
addpath(sprintf('%s/cbrewer/',lib_path));
addpath(sprintf('%s/export_fig/',lib_path));
addpath(sprintf('%s/ATSA_Meko/',lib_path));
addpath(sprintf('%s/freezeColors/',lib_path));
defcolormap;

% --------------------------
%         Setting  
% --------------------------
global domain_def;
global chosen_prd;
global nyr nyr_gs;
global nds ds_s ds_e nmd mdi;
global stime;
global TT;
global DATA_VI_05rs_out DATA_VI_05rs_sm_out;
global DATA_VI_05rs_bw_out DATA_VI_05rs_bw_sm_out;
global DATA_VI_05rs_out_gs;
global DATA_VI_05rs_out_gs1;
global DATA_VI_05rs_bw_out;
global DATA_VI_05rs_dtds_out;
global DATA_VI_05rs_dtds_out_gs;
global DATA_VI_05rs_dtds_out_gs1;
global DATA_VI_05rs_dtdspw_out;
global mon_ann_r mon_ann_p;
global mon_ann_r_pw mon_ann_p_pw;
global mw_r mw_p;
global mw_r_pw mw_p_pw;
global r_trnd r_trnd_p;
global r_trnd_pw;
global mon_anm_ratio mon_anm_ratio_mean;
global mw_ratio;
global ratio_trnd ratio_trnd_p;
global mon_cv mon_std;
global mw_cv mw_std;
global cv_trnd cv_trnd_p std_trnd std_trnd_p;

global window_wd nw nw_gs;

global DATA_ERA5_05rs_out;
global DATA_CRU_05rs_out;
global sos_map eos_map lgs_map;

global DATA_KG_CLMZ_05rs_out;
global DATA_LC_out;
global lc_dom_grp clmzone_grp;
global DATA_SMRoot_05rs_out;

global DATA_CLM_05rs_out_gs DATA_CLM_05rs_out_gs1;
global DATA_CLM_05rs_dtds_out DATA_CLM_05rs_dtds_out_gs1;
global DATA_CLM_05rs_dtdspw_out;

global mon_vegclm_par_r mon_vegclm_par_p;
global DATA_CLMINX_out DATA_CLMINX_out_gs DATA_CLMINX_out_gs1;
global v_clminx_r v_clminx_p;

global DATA_CLM_05rs_gs_lag_out DATA_CLM_05rs_gs1_lag_out;
global DATA_CLM_05rs_dtds_out_gs1;
global ax_show;
global clm_cl_x;

% Domain definition, generate data mask
%            Lon.W., Lon.E., Lat.S, Lat.N
domain_def = [-13.0 33.0 34.0 72.5 ];
chosen_mons1=[4:8]; % Growing season mean vegetation
chosen_prd=1982:2014;
nyr=chosen_prd(end)-chosen_prd(1)+1;
nyr_gs=nyr-1;

% Monthly calendar
stime = datetime(chosen_prd(1),1,15);
tstep = calmonths(1);
TT = timetable('Size',[nyr*12 1],'VariableTypes',{'double'},...
               'TimeStep',tstep,'StartTime',stime);
ds_s=1; ds_e=2; nds=ds_e-ds_s+1; 
nmd=4;
mdi=4; % Double logistic from TIMESAT

% Parameters for calculating trend
global window_wd nw;
window_wd=20;
nw=(nyr-window_wd)+1;
nw_gs=nw-1;

% Setting for plotting
ax_show=struct('south_s',domain_def(3),'north_s',domain_def(4),'west_s',domain_def(1),'east_s',domain_def(2),'resolution_s',0.5, ...
               'south_d',domain_def(3),'north_d',domain_def(4),'west_d',domain_def(1),'east_d',domain_def(2),'resolution_d',0.5, ...
               'countries',0,'contourf', 0, 'lake',1,'coast',1,'lonlat',0,'frame',1, 'landmask',0);

global lons lats;
[lons,lats] = meshgrid(ax_show.south_s:0.5:ax_show.north_s, ax_show.west_s:0.5:ax_show.east_s);
lons=lons(1:(end-1),1:(end-1));
lats=lats(1:(end-1),1:(end-1));
%% ------------------------------
%           Climate
% ------------------------------
% ---------- input -------------
 preproc_ERA5; % <--
 preproc_CRU;  % <--
 preproc_climatezone; % <--
 preproc_MODIS_Landcover_Global; % <--
 preproc_GLEAM;  % <--
%% ------------------------------
%          Vegetation 
% ------------------------------
% ---------- input -------------
 preproc_GIMMS_NDVI3g_Global('RS_05'); % <-- Paper I: NDVI
 preproc_VIP_NDVI_EVI2_Global('RS_05'); % <-- Paper I: EVI2
%  preproc_MODIS_NDVI_EVI_Global('RS_05'); % <-- Paper I: NDSOS_EOS_LGS() % <-- Paper I >>> (Fig. S1)VI, EVI
 % preproc_SPOT_LAI_Global('RS_05'); % <-- Paper I: LAI
 % preproc_FluxCom_Global('RS_05'); % <-- Paper I
 
% Import TIMESAT smoothed data
 TIMESAT_smoothed;
  
 % NDVI, EVI: positive values generally indicate the presence of vegetation 
 % (with greater values indicating healthier vegetation). 
 % Negative values generally indicate a lack of vegetation (water, rock, soil, etc.).
 % LAI should > 0
 DATA_VI_05rs_out(DATA_VI_05rs_out<=0)=nan;
 DATA_VI_05rs_sm_out(DATA_VI_05rs_sm_out<=0)=nan;
 
% use mean to make grid point nan if any of the grid point 
% (spatially and temporally) with nan occurs
global mask_gs;
mask_gs=squeeze(mean(mean(mean(DATA_VI_05rs_sm_out,3),4),5));
mask_gs(~isnan(mask_gs))=1;

% Define growing season
SOS_EOS_LGS() % <-- Paper I >>> (Fig. S1)

% ---------- calculate ------------- 
preproc_calcVI_I('monthly_dsclim',0,1); % <-- Paper I:  % 1st param: 1: prewhitened, 0: orig <-- Paper I
                                                        % 2nd param: 1: use growing season definition, 0: do not use
%%
month_to_tot_month_anom_calc(0,1); %  <-- Paper I: % 1st param: 1: prewhitened, 0: orig
                                                   % 2nd param: 1: use growing season definition, 0: do not use
month_to_tot_month_anom_calc_ratio(1);  % <-- Paper I: 1: use growing season definition, 0: do not use
month_to_tot_month_anom_calc_cv(1);  % <-- Paper I: 1: use growing season definition, 0: do not use

% loadsave('save','WORKSAPCE');

% ---------- viusalize -------------
% month_to_tot_month_anom_vsl;
% max_ssn_iav; % <-- Paper I >>> (Fig. 1)
% max_ssn_iav_ratio; % <-- Paper I >>> (Fig. Sx)

max_ssn_iav_anom; % <-- Paper I >>> (Fig. 1)

%% Land class and climate zone
iav_by_landclass; % <-- Paper I >>> (Fig. Sx)

% Stability map and rank by land class, climate zone
cv_by_lc_clmzone(); % <-- Paper I >>> (Fig. 2)

%% ------------------------------
%          Local climate
% ------------------------------
preproc_calcLocalClim('monthly_dsclim',0,1);

decomp_VegClm_relation();
Maxwell_color_triangle(31);
decomp_VegClm_relation_plot();

%%
%         Teleconnection
% ------------------------------
preproc_ClmInx(1);
var_to_clminx_month_calc_r(); % Climate teleconnection
dom_clminx_plot();
dom_climinx_supp_plot();

% Analysis_framework();
% smoothing_methods();

%%
%    Vegetation and lagged local climate
% ---------------------------------------




% --- Decompose variability
% EOF;
% SVD; % <--

%% Load or Save data
% DATA_CALC_VI_05RS_pw DATA_CALC_VI_05RS, FLUXCOM, VIP_EVI2, VIP_NDVI, NDVI3g, MODIS, CLIMATE, WORKSAPCE
% loadsave('save','DATA_CALC_VI_05RS');

%% Idea of increasing control of spring anomalies
% Calculated data from preproc_calcVI_I.m
% act_type='load';
% save_type={'NDVI3g'}; % NDVI3g, MODIS, CLIMATE, WORKSAPCE

% Calculated data from spring_to_tot_month_anom.m
% load(sprintf('%s/Calc_sp_var_ann.mat',wrk_dir));

% visualize 
% spring_to_tot_month_anom_vsl



