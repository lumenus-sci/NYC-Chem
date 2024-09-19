seqs = [x for x in range(1,49)]
sdays = [*(x for x in range(20,32) for _ in (0,1)), *(x for x in range(1,13) for _ in (0,1))]
edays = [20, *(x for x in range(21,32) for _ in (0,1)), *(x for x in range(1,13) for _ in (0,1)), 13]
smons = [*([7] * 24), *([8] * 24)]
emons = [*([7] * 23), *([8] * 25)]
shrs = [0, 12] * 24
ehrs = [12, 0] * 24


for seq, smon, sday, shr, emon, eday, ehr in zip(seqs, smons, sdays, shrs, emons, edays, ehrs):
    with open(f'wrf_ghg_{seq:02}.sh','a') as fl:
        fl.write(f'&time_control\n')
        fl.write(f'run_days                            = 0,\n')
        fl.write(f'run_hours                           = 12,\n')
        fl.write(f'run_minutes                         = 0,\n')
        fl.write(f'run_seconds                         = 0,\n')
        fl.write(f'start_year                          = 2023, 2023, 2023, 2023,\n')
        fl.write(f'start_month                         = {smon:02}, {smon:02}, {smon:02}, {smon:02},\n')
        fl.write(f'start_day                           = {sday:02}, {sday:02}, {sday:02}, {sday:02},\n')
        fl.write(f'start_hour                          = {shr:02}, {shr:02}, {shr:02}, {shr:02},\n')
        fl.write(f'start_minute                        = 00, 00, 00, 00,\n')
        fl.write(f'start_second                        = 00, 00, 00, 00,\n')
        fl.write(f'end_year                            = 2023, 2023, 2023, 2023,\n')
        fl.write(f'end_month                           = {emon:02}, {emon:02}, {emon:02}, {emon:02},\n')
        fl.write(f'end_day                             = {eday:02}, {eday:02}, {eday:02}, {eday:02},\n')
        fl.write(f'end_hour                            = {ehr:02}, {ehr:02}, {ehr:02}, {ehr:02},\n')
        fl.write(f'end_minute                          = 00, 00, 00, 00,\n')
        fl.write(f'end_second                          = 00, 00, 00, 00,\n')
        fl.write(f'interval_seconds                    = 21600\n')
        fl.write(f'auxinput4_interval                  = 360,360,360,360,360\n')
        fl.write(f'auxinput4_inname                    = "wrflowinp_d<domain>"\n')
        fl.write(f'input_from_file                     = .true.,.true.,.true.,.true.,\n')
        fl.write(f'history_interval                    = 60, 60, 60, 60,\n')
        fl.write(f'frames_per_outfile                  = 1,  1,  1,  1,\n')
        if seq == 1:
            fl.write(f'restart                             = .false.,\n')
        else:
            fl.write(f'restart                             = .true.,\n')
        fl.write(f'restart_interval                    = 720,\n')
        fl.write(f'auxinput5_interval_d                = 31, 31, 31, 31, 31,\n')
        fl.write(f'auxinput15_interval_d               = 8, 8, 8, 8, 8,\n')
        fl.write(f'io_form_history                     = 2\n')
        fl.write(f'io_form_restart                     = 2\n')
        fl.write(f'io_form_input                       = 2\n')
        fl.write(f'io_form_boundary                    = 2\n')
        fl.write(f'io_form_auxinput4                   = 2\n')
        fl.write(f'io_form_auxinput5                   = 2\n')
        fl.write(f'io_form_auxinput6                   = 2\n')
        fl.write(f'io_form_auxinput12                  = 2\n')
        fl.write(f'io_form_auxinput15                  = 2\n')
        fl.write(f"auxinput15_inname                    = 'wrfVPRM_parameterFromMODIS_d<domain>',\n")
        fl.write(f'force_use_old_data                  = .true.\n')
        fl.write(f"auxinput12_inname                    = 'wrf_chem_input',\n")
        fl.write(f"iofields_filename                   = 'hist_io_mods_d01', 'hist_io_mods_d01', 'hist_io_mods_d01', 'hist_io_mods_d01',\n")
        fl.write(f'debug_level                         = 0\n')
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&domains\n')
        fl.write(f'time_step                           = 60,\n')
        fl.write(f'time_step_fract_num                 = 0,\n')
        fl.write(f'time_step_fract_den                 = 1,\n')
        fl.write(f'max_dom                             = 2,\n')
        fl.write(f's_we                                = 1,     1,     1,     1,\n')
        fl.write(f'e_we                                = 443,   301,   161,   161,\n')
        fl.write(f's_sn                                = 1,     1,     1,     1,\n')
        fl.write(f'e_sn                                = 266,   256,   161,   161,\n')
        fl.write(f's_vert                              = 1,     1,     1,     1,\n')
        fl.write(f'e_vert                              = 48,   48,   48,     48,\n')
        fl.write(f'sfcp_to_sfcp                         = .false.\n')
        fl.write(f'dx                                  = 12000,   800,   160,   32,\n')
        fl.write(f'dy                                  = 12000,   800,   160,   32,\n')
        fl.write(f'grid_id                             = 1,     2,     3,     4,\n')
        fl.write(f'parent_id                           = 1,   1,   2,   3,\n')
        fl.write(f'i_parent_start                      = 1,   372,   164,   77,\n')
        fl.write(f'j_parent_start                      = 1,   154,   90,   80,\n')
        fl.write(f'parent_grid_ratio                   = 1,   15,   5,   5,\n')
        fl.write(f'parent_time_step_ratio              = 1,   15,    5,    5,\n')
        fl.write(f'feedback                            = 0,\n')
        fl.write(f'smooth_option                       = 0,\n')
        fl.write(f'num_metgrid_levels                  = 34,\n')
        fl.write(f'num_metgrid_soil_levels             = 4,\n')
        fl.write(f'interp_type                         = 2,\n')
        fl.write(f'lagrange_order                      = 2,\n')
        fl.write(f'zap_close_levels                    = 500,\n')
        fl.write(f'lowest_lev_from_sfc                 = .false.,\n')
        fl.write(f'force_sfc_in_vinterp                = 1,\n')
        fl.write(f'p_top_requested                     = 1000,\n')
        fl.write(f'eta_levels                          = 1.000,0.997,0.994,0.991,0.988,0.985,0.980,0.975,0.970,0.960,0.950,\n') 
        fl.write(f'                                    0.940,0.930,0.920,0.910,0.895,0.880,\n')
        fl.write(f'                                    0.865,0.850,0.825,0.800,0.775,0.750,\n')
        fl.write(f'                                    0.720,0.690,0.660,0.630,0.600,0.570,\n')
        fl.write(f'                                    0.540,0.510,0.475,0.440,0.405,0.370,\n')
        fl.write(f'                                    0.330,0.290,0.250,0.210,0.175,0.145,\n')
        fl.write(f'                                    0.115,0.090,0.065,0.045,0.025,0.010,\n')
        fl.write(f'                                    0.000,\n')
        fl.write(f'smooth_option                       = 0\n')
        fl.write(f'vortex_interval                     = 15\n')
        fl.write(f'max_vortex_speed                    = 80\n')
        fl.write(f'corral_dist                         = 125\n')
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&physics\n')
        fl.write(f'mp_physics                          = 10,     10,     10,     10,     10,\n')
        fl.write(f'do_radar_ref           = 1 ,\n')
        fl.write(f'ra_lw_physics                       = 1,     1,     1,     1,\n') 
        fl.write(f'ra_sw_physics                       = 1,     1,     1,     1,\n') 
        fl.write(f'radt                                = 30,    30,    30,   30,\n') 
        fl.write(f'sf_sfclay_physics                   = 1,     1,     1,     1,\n')
        fl.write(f'sf_surface_physics                  = 2,     2,     2,     2,\n')
        fl.write(f'bl_pbl_physics                      = 11,    11,     0,     0,\n')
        fl.write(f'bldt                                = 0,     0,     0,     0,\n')
        fl.write(f'cu_physics                          = 3,     3,     0,     0,\n')
        fl.write(f'cu_rad_feedback                     = .true.\n')
        fl.write(f'cu_diag                             = 1,     1,     0,     0,\n')
        fl.write(f'cudt                                = 5,     5,     5,     5,\n')
        fl.write(f'isfflx                              = 1,\n')
        fl.write(f'ifsnow                              = 0,\n')
        fl.write(f'icloud                              = 1,\n')
        fl.write(f'surface_input_source                = 1,\n')
        fl.write(f'sst_update                          = 0, ! if input_from_file false, sst_update needs to be 0\n')
        fl.write(f'num_soil_layers                     = 4,\n')
        fl.write(f'maxiens                             = 1,\n')
        fl.write(f'maxens                              = 3,\n')
        fl.write(f'maxens2                             = 3,\n')
        fl.write(f'maxens3                             = 16,\n')
        fl.write(f'ensdim                              = 144,\n')
        fl.write(f'prec_acc_dt                         = 1440,  1440,  1440,  1440,\n') 
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&fdda\n')
        fl.write(f'grid_fdda                           = 1,       1,      1,      1,\n')
        fl.write(f'gfdda_inname                        = "wrffdda_d<domain>",\n')
        fl.write(f'gfdda_end_h                         = 12,      12,     12,     12,\n')
        fl.write(f'gfdda_interval_m                    = 360,     720,    360,    360,\n')
        fl.write(f'grid_sfdda                          = 0,       0,      0,      0,\n')
        fl.write(f'sgfdda_inname                       = "wrfsfdda_d<domain>",\n')
        fl.write(f'sgfdda_end_h                        = 96,      96,     96,     96,\n')
        fl.write(f'sgfdda_interval_m                   = 360,     360,    360,    360,\n')
        fl.write(f'rinblw                              = 250,\n')
        fl.write(f'guv_sfc                             = 0.0003, 0.0001,\n')
        fl.write(f'fgdt                                = 0,       0,      0,      0,\n')
        fl.write(f'if_no_pbl_nudging_uv                = 1,       1,      1,      1,\n')
        fl.write(f'if_no_pbl_nudging_t                 = 1,       1,      1,      1,\n')
        fl.write(f'if_no_pbl_nudging_q                 = 1,       1,      1,      1,\n')
        fl.write(f'if_no_pbl_nudging_ph                = 1,       1,      1,      1,\n')
        fl.write(f'if_zfac_uv                          = 0,       0,      0,      0,\n')
        fl.write(f'k_zfac_uv                          = 10,      10,     1,      1,\n')
        fl.write(f'if_zfac_t                           = 0,       0,      0,      0,\n')
        fl.write(f'k_zfac_t                           = 10,      10,     1,      1,\n')
        fl.write(f'if_zfac_q                           = 0,       0,      0,      0,\n')
        fl.write(f'k_zfac_q                           = 10,      10,     1,      1,\n')
        fl.write(f'if_zfac_ph                          = 0,       0,      0,      0,\n')
        fl.write(f'k_zfac_ph                           = 10,      10,      10,      10,\n')
        fl.write(f'guv                                 = 3.0E-5,  3.0E-5, 3.0E-5, 3.0E-5,\n')
        fl.write(f'gph                                 = 3.0E-5,  3.0E-5, 3.0E-5, 3.0E-5,\n')
        fl.write(f'gt                                  = 3.0E-5,  3.0E-5, 3.0E-5, 3.0E-5,\n')
        fl.write(f'gq                                  = 3.0E-5,  1.0E-5, 0,      0,\n')
        fl.write(f'if_ramping                          = 1,\n')
        fl.write(f'xwavenum                            = 5,   1,\n')  
        fl.write(f'ywavenum                            = 3,   1,\n') 
        fl.write(f'dtramp_min                          = 60.0,\n')
        fl.write(f'io_form_gfdda                       = 2,\n')
        fl.write(f'io_form_sgfdda                      = 2,\n')
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&dynamics\n')
        fl.write(f'w_damping                           = 1,\n')
        fl.write(f'diff_opt                            = 0, 0, 2, 2,\n')
        fl.write(f'km_opt                              = 1, 1, 2, 2,\n')
        fl.write(f'diff_6th_opt                        = 0,\n')
        fl.write(f'diff_6th_factor                     = 0.12,\n')
        fl.write(f'damp_opt                            = 3,\n')
        fl.write(f'base_temp                           = 290.\n')
        fl.write(f'zdamp                               = 5000.,  5000.,  5000.,  5000.,\n')
        fl.write(f'dampcoef                            = 0.2,    0.2,    0.2,    0.2,\n')
        fl.write(f'khdif                               = 0,      0,      0,      0,\n')
        fl.write(f'kvdif                               = 0,      0,      0,      0,\n')
        fl.write(f'smdiv                               = 0.1,    0.1,    0.1,    0.1,\n')
        fl.write(f'emdiv                               = 0.01,   0.01,   0.01,   0.01,\n')
        fl.write(f'epssm                               = 0.1,    0.1,    0.1     0.1,\n')
        fl.write(f'non_hydrostatic                     = .true., .true., .true., .true.,\n')
        fl.write(f'time_step_sound                     = 4,      4,      4,      4,\n')
        fl.write(f'chem_adv_opt                        = 2,      2,       2,       2,       2,\n')
        fl.write(f'moist_adv_opt                       = 2,      2,       2,       2,       2,\n')
        fl.write(f'scalar_adv_opt                      = 2,      2,       2,       2,       2,\n')
        fl.write(f'tke_adv_opt                         = 2,      2,       2,       2,       2,\n')
        fl.write(f'tracer_adv_opt                      = 2,      2,      2,       2,       2,\n')
        fl.write(f'h_mom_adv_order                     = 5,      5,      5,      5,\n')
        fl.write(f'v_mom_adv_order                     = 3,      3,      3,      3,\n')
        fl.write(f'h_sca_adv_order                     = 5,      5,      5,      5,\n')
        fl.write(f'v_sca_adv_order                     = 3,      3,      3,      3,\n')
        fl.write(f'do_avgflx_em                        = 1, 1, 1, 1, 1, 1, 1, 1, 1,\n')
        fl.write(f'do_avgflx_cugd                      = 1, 1, 1, 1, 1, 1, 1, 1, 1,\n')
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&bdy_control\n')
        fl.write(f'spec_bdy_width                      = 5,\n')
        fl.write(f'spec_zone                           = 1,\n')
        fl.write(f'relax_zone                          = 4,\n')
        fl.write(f'specified                           = .true., .false.,.false.,.false.,\n')
        fl.write(f'periodic_x                          = .false.,.false.,.false.,.false.,\n')
        fl.write(f'symmetric_xs                        = .false.,.false.,.false.,.false.,\n')
        fl.write(f'symmetric_xe                        = .false.,.false.,.false.,.false.,\n')
        fl.write(f'open_xs                             = .false.,.false.,.false.,.false.,\n')
        fl.write(f'open_xe                             = .false.,.false.,.false.,.false.,\n')
        fl.write(f'periodic_y                          = .false.,.false.,.false.,.false.,\n')
        fl.write(f'symmetric_ys                        = .false.,.false.,.false.,.false.,\n')
        fl.write(f'symmetric_ye                        = .false.,.false.,.false.,.false.,\n')
        fl.write(f'open_ys                             = .false.,.false.,.false.,.false.,\n')
        fl.write(f'open_ye                             = .false.,.false.,.false.,.false.,\n')
        fl.write(f'nested                              = .false., .true., .true., .true.,\n')
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&grib2\n')
        fl.write(f'/\n')
        fl.write(f'&chem\n')
        fl.write(f'kemit                               = 12,\n')
        fl.write(f'chem_opt                            = 17,      17,      17,      17,      17,\n')
        fl.write(f"vprm_opt                            = 'VPRM_table_US', 'VPRM_table_US', 'VPRM_table_US','VPRM_table_US','VPRM_table_US',\n")
        fl.write(f'bioemdt                             = 30,       30,       30,       30,       30,\n')
        fl.write(f'photdt                              = 30,       30,       30,       30,       30,\n')
        fl.write(f'chemdt                              = 1.,       1.,       1.,       1.,       1.,\n')
        fl.write(f'io_style_emissions                  = 1,\n')
        fl.write(f'emiss_inpt_opt                      = 16,        16,        16,        16,        16,\n')
        fl.write(f'emiss_opt                           = 17,        17,        17,        17,        17,\n')
        fl.write(f'chem_in_opt                         = 0,        0,        0,        0,        0,\n')
        fl.write(f'phot_opt                            = 1,        1,        1,       1,       1,\n')
        fl.write(f'gas_drydep_opt                      = 0,        0,        0,        0,        0,\n')
        fl.write(f'bio_emiss_opt                       = 17,        17,        17,        17,        17,\n')
        fl.write(f'ne_area                             = 60,\n')
        fl.write(f'dust_opt                            = 0,\n')
        fl.write(f'dmsemis_opt                         = 0,\n')
        fl.write(f'seas_opt                            = 0,\n')
        fl.write(f'gas_bc_opt                          = 1,        1,        1,        1,        1,\n')
        fl.write(f'gas_ic_opt                          = 1,        1,        1,        1,        1,\n')
        fl.write(f'aer_bc_opt                          = 1,        1,        1,        1,        1,\n')
        fl.write(f'aer_ic_opt                          = 1,        1,        1,        1,        1,\n')
        fl.write(f'gaschem_onoff                       = 1,        1,        1,        1,        1,\n')
        fl.write(f'aerchem_onoff                       = 1,        1,        1,        1,        1,\n')
        fl.write(f'wetscav_onoff                       = 0,        0,        0,        0,        0,\n')
        fl.write(f'cldchem_onoff                       = 0,        0,        0,        0,        0,\n')
        fl.write(f'vertmix_onoff                       = 1,        1,        1,        1,        1,\n')
        fl.write(f'chem_conv_tr                        = 1,        0,        0,        0,        0,\n')
        fl.write(f'conv_tr_wetscav                     = 0, \n')
        fl.write(f'biomass_burn_opt                    = 0,        0,        0,        0,        0,\n')
        fl.write(f'plumerisefire_frq                   = 30,       30,       30,       30,       30,\n')
        fl.write(f'aer_ra_feedback                     = 0,        0,        0,        0,        0,\n')
        fl.write(f'have_bcs_chem                       = .true., .true., .true., .true., .true.,\n')
        fl.write(f'have_bcs_tracer                     = .true., .true., .true., .true., .true.,\n')
        fl.write(f'/\n')
        fl.write(f'\n')
        fl.write(f'&namelist_quilt\n')
        fl.write(f'nio_tasks_per_group = 0,\n')
        fl.write(f'nio_groups = 2,\n')
        fl.write(f'/\n')

