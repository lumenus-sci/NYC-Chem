seqs = [x for x in range(1,49)]
days = [x for x in range(20,32) for _ in (0,1)]
days2 = [x for x in range(1,13) for _ in (0,1)]
days.extend(days2)
mons = [7] * 24
mons2 = [8] * 24
mons.extend(mons2)
hrs = [0, 12] * 24


for i, seq, mon, day, hr in enumerate(zip(seqs, mons, days, hrs)):
    with open(f'wrf_ghg_{seq:02}.sh','a') as fl:
        &time_control\n')
        run_days                            = 0,\n')
        run_hours                           = 12,\n')
        run_minutes                         = 0,\n')
        run_seconds                         = 0,\n')
        start_year                          = 2023, 2023, 2023, 2023,\n')
        start_month                         = 07, 07, 07, 07,\n')
        start_day                           = 20, 20, 20, 20,\n')
        start_hour                          = 00, 00, 00, 00,\n')
        start_minute                        = 00, 00, 00, 00,\n')
        start_second                        = 00, 00, 00, 00,\n')
        end_year                            = 2023, 2023, 2023, 2023,\n')
        end_month                           = 07, 07, 07, 07,\n')
        end_day                             = 20, 20, 20, 20,\n')
        end_hour                            = 12, 12, 12, 12,\n')
        end_minute                          = 00, 00, 00, 00,\n')
        end_second                          = 00, 00, 00, 00,\n')
        interval_seconds                    = 21600\n')
        auxinput4_interval                  = 360,360,360,360,360\n')
        auxinput4_inname                    = "wrflowinp_d<domain>"\n')
        input_from_file                     = .true.,.true.,.true.,.true.,\n')
        history_interval                    = 60, 60, 60, 60,\n')
        frames_per_outfile                  = 1,  1,  1,  1,\n')
        restart                             = .false.,\n')
        restart_interval                    = 720,\n')
        auxinput5_interval_d                = 31, 31, 31, 31, 31,\n')
        auxinput15_interval_d               = 8, 8, 8, 8, 8,\n')
        io_form_history                     = 2\n')
        io_form_restart                     = 2\n')
        io_form_input                       = 2\n')
        io_form_boundary                    = 2\n')
        io_form_auxinput4                   = 2\n')
        io_form_auxinput5                   = 2\n')
        io_form_auxinput6                   = 2\n')
        io_form_auxinput12                  = 2\n')
        io_form_auxinput15                  = 2\n')
        auxinput15_inname                    = 'wrfVPRM_parameterFromMODIS_d<domain>',\n')
        force_use_old_data                  = .true.\n')
        auxinput12_inname                    = 'wrf_chem_input',\n')
        iofields_filename                   = 'hist_io_mods_d01', 'hist_io_mods_d01', 'hist_io_mods_d01', 'hist_io_mods_d01',\n')
        debug_level                         = 0\n')
        /\n')
        \n')
        &domains\n')
        time_step                           = 60,\n')
        time_step_fract_num                 = 0,\n')
        time_step_fract_den                 = 1,\n')
        max_dom                             = 2,\n')
        s_we                                = 1,     1,     1,     1,\n')
        e_we                                = 443,   301,   161,   161,\n')
        s_sn                                = 1,     1,     1,     1,\n')
        e_sn                                = 266,   256,   161,   161,\n')
        s_vert                              = 1,     1,     1,     1,\n')
        e_vert                              = 48,   48,   48,     48,\n')
        sfcp_to_sfcp                         = .false.\n')
        dx                                  = 12000,   800,   160,   32,\n')
        dy                                  = 12000,   800,   160,   32,\n')
        grid_id                             = 1,     2,     3,     4,\n')
        parent_id                           = 1,   1,   2,   3,\n')
        i_parent_start                      = 1,   372,   164,   77,\n')
        j_parent_start                      = 1,   154,   90,   80,\n')
        parent_grid_ratio                   = 1,   15,   5,   5,\n')
        parent_time_step_ratio              = 1,   15,    5,    5,\n')
        feedback                            = 0,\n')
        smooth_option                       = 0,\n')
        num_metgrid_levels                  = 34,\n')
        num_metgrid_soil_levels             = 4,\n')
        interp_type                         = 2,\n')
        lagrange_order                      = 2,\n')
        zap_close_levels                    = 500,\n')
        lowest_lev_from_sfc                 = .false.,\n')
        force_sfc_in_vinterp                = 1,\n')
        p_top_requested                     = 1000,\n')
        eta_levels                          = 1.000,0.997,0.994,0.991,0.988,0.985,0.980,0.975,0.970,0.960,0.950,\n') 
                                            0.940,0.930,0.920,0.910,0.895,0.880,\n')
                                            0.865,0.850,0.825,0.800,0.775,0.750,\n')
                                            0.720,0.690,0.660,0.630,0.600,0.570,\n')
                                            0.540,0.510,0.475,0.440,0.405,0.370,\n')
                                            0.330,0.290,0.250,0.210,0.175,0.145,\n')
                                            0.115,0.090,0.065,0.045,0.025,0.010,\n')
                                            0.000,\n')
        smooth_option                       = 0\n')
        vortex_interval                     = 15\n')
        max_vortex_speed                    = 80\n')
        corral_dist                         = 125\n')
        /\n')
        \n')
        &physics\n')
        mp_physics                          = 10,     10,     10,     10,     10,\n')
        do_radar_ref           = 1 ,\n')
        ra_lw_physics                       = 1,     1,     1,     1,\n') 
        ra_sw_physics                       = 1,     1,     1,     1,\n') 
        radt                                = 30,    30,    30,   30,\n') 
        sf_sfclay_physics                   = 1,     1,     1,     1,\n')
        sf_surface_physics                  = 2,     2,     2,     2,\n')
        bl_pbl_physics                      = 11,    11,     0,     0,\n')
        bldt                                = 0,     0,     0,     0,\n')
        cu_physics                          = 3,     3,     0,     0, 
        cu_rad_feedback                     = .true.
        cu_diag                             = 1,     1,     0,     0,
        cudt                                = 5,     5,     5,     5, 
        isfflx                              = 1,
        ifsnow                              = 0,
        icloud                              = 1,
        surface_input_source                = 1,
        sst_update                          = 0, ! if input_from_file false, sst_update needs to be 0
        num_soil_layers                     = 4,
        maxiens                             = 1,
        maxens                              = 3,
        maxens2                             = 3,
        maxens3                             = 16,
        ensdim                              = 144,
        prec_acc_dt                         = 1440,  1440,  1440,  1440, 
        /

        &fdda
        grid_fdda                           = 0,       0,      0,      0,
        gfdda_inname                        = "wrffdda_d<domain>",
        gfdda_end_h                         = 12,      12,     12,     12,
        gfdda_interval_m                    = 360,     720,    360,    360,
        grid_sfdda                          = 0,       0,      0,      0,
        sgfdda_inname                       = "wrfsfdda_d<domain>",
        sgfdda_end_h                        = 96,      96,     96,     96, 
        sgfdda_interval_m                   = 360,     360,    360,    360,
        rinblw                              = 250,                     
        guv_sfc                             = 0.0003, 0.0001,
        fgdt                                = 0,       0,      0,      0,
        if_no_pbl_nudging_uv                = 1,       1,      1,      1,
        if_no_pbl_nudging_t                 = 1,       1,      1,      1,
        if_no_pbl_nudging_q                 = 1,       1,      1,      1,
        if_no_pbl_nudging_ph                = 1,       1,      1,      1, 
        if_zfac_uv                          = 0,       0,      0,      0,
        k_zfac_uv                          = 10,      10,     1,      1,
        if_zfac_t                           = 0,       0,      0,      0,
        k_zfac_t                           = 10,      10,     1,      1,
        if_zfac_q                           = 0,       0,      0,      0,
        k_zfac_q                           = 10,      10,     1,      1,
        if_zfac_ph                          = 0,       0,      0,      0,
        k_zfac_ph                           = 10,      10,      10,      10,
        guv                                 = 3.0E-5,  3.0E-5, 3.0E-5, 3.0E-5,
        gph                                 = 3.0E-5,  3.0E-5, 3.0E-5, 3.0E-5, 
        gt                                  = 3.0E-5,  3.0E-5, 3.0E-5, 3.0E-5, 
        gq                                  = 3.0E-5,  1.0E-5, 0,      0,
        if_ramping                          = 1,
        xwavenum                            = 5,   1,  
        ywavenum                            = 3,   1,  
        dtramp_min                          = 60.0,
        io_form_gfdda                       = 2,
        io_form_sgfdda                      = 2,
        /


        &dynamics
        w_damping                           = 1,
        diff_opt                            = 0, 0, 2, 2,
        km_opt                              = 1, 1, 2, 2,
        diff_6th_opt                        = 0,
        diff_6th_factor                     = 0.12,
        damp_opt                            = 3,
        base_temp                           = 290.
        zdamp                               = 5000.,  5000.,  5000.,  5000.,
        dampcoef                            = 0.2,    0.2,    0.2,    0.2,
        khdif                               = 0,      0,      0,      0,
        kvdif                               = 0,      0,      0,      0,
        smdiv                               = 0.1,    0.1,    0.1,    0.1,
        emdiv                               = 0.01,   0.01,   0.01,   0.01,
        epssm                               = 0.1,    0.1,    0.1     0.1,
        non_hydrostatic                     = .true., .true., .true., .true.,
        time_step_sound                     = 4,      4,      4,      4,
        chem_adv_opt                        = 2,      2,       2,       2,       2,
        moist_adv_opt                       = 2,      2,       2,       2,       2,
        scalar_adv_opt                      = 2,      2,       2,       2,       2,
        tke_adv_opt                         = 2,      2,       2,       2,       2,
        tracer_adv_opt                      = 2,      2,      2,       2,       2,
        h_mom_adv_order                     = 5,      5,      5,      5,
        v_mom_adv_order                     = 3,      3,      3,      3,
        h_sca_adv_order                     = 5,      5,      5,      5,
        v_sca_adv_order                     = 3,      3,      3,      3,
        do_avgflx_em                        = 1, 1, 1, 1, 1, 1, 1, 1, 1,
        do_avgflx_cugd                      = 1, 1, 1, 1, 1, 1, 1, 1, 1,
        /

        &bdy_control
        spec_bdy_width                      = 5,
        spec_zone                           = 1,
        relax_zone                          = 4,
        specified                           = .true., .false.,.false.,.false.,
        periodic_x                          = .false.,.false.,.false.,.false.,
        symmetric_xs                        = .false.,.false.,.false.,.false.,
        symmetric_xe                        = .false.,.false.,.false.,.false.,
        open_xs                             = .false.,.false.,.false.,.false.,
        open_xe                             = .false.,.false.,.false.,.false.,
        periodic_y                          = .false.,.false.,.false.,.false.,
        symmetric_ys                        = .false.,.false.,.false.,.false.,
        symmetric_ye                        = .false.,.false.,.false.,.false.,
        open_ys                             = .false.,.false.,.false.,.false.,
        open_ye                             = .false.,.false.,.false.,.false.,
        nested                              = .false., .true., .true., .true.,
        /

        &grib2
        /
        &chem
        kemit                               = 12,
        chem_opt                            = 17,      17,      17,      17,      17,
        vprm_opt                            = 'VPRM_table_US', 'VPRM_table_US', 'VPRM_table_US','VPRM_table_US','VPRM_table_US',
        bioemdt                             = 30,       30,       30,       30,       30,
        photdt                              = 30,       30,       30,       30,       30,
        chemdt                              = 1.,       1.,       1.,       1.,       1.,
        io_style_emissions                  = 1,
        emiss_inpt_opt                      = 16,        16,        16,        16,        16,
        emiss_opt                           = 17,        17,        17,        17,        17,
        chem_in_opt                         = 0,        0,        0,        0,        0,
        phot_opt                            = 1,        1,        1,       1,       1,
        gas_drydep_opt                      = 0,        0,        0,        0,        0,
        bio_emiss_opt                       = 17,        17,        17,        17,        17,
        ne_area                             = 60,
        dust_opt                            = 0,
        dmsemis_opt                         = 0,
        seas_opt                            = 0,
        gas_bc_opt                          = 1,        1,        1,        1,        1,
        gas_ic_opt                          = 1,        1,        1,        1,        1,
        aer_bc_opt                          = 1,        1,        1,        1,        1,
        aer_ic_opt                          = 1,        1,        1,        1,        1,
        gaschem_onoff                       = 1,        1,        1,        1,        1,
        aerchem_onoff                       = 1,        1,        1,        1,        1,
        wetscav_onoff                       = 0,        0,        0,        0,        0,
        cldchem_onoff                       = 0,        0,        0,        0,        0,
        vertmix_onoff                       = 1,        1,        1,        1,        1,
        chem_conv_tr                        = 1,        0,        0,        0,        0,
        conv_tr_wetscav                     = 0, 
        biomass_burn_opt                    = 0,        0,        0,        0,        0,
        plumerisefire_frq                   = 30,       30,       30,       30,       30,
        aer_ra_feedback                     = 0,        0,        0,        0,        0,
        have_bcs_chem                       = .true., .true., .true., .true., .true.,
        have_bcs_tracer                     = .true., .true., .true., .true., .true.,
        /

        &namelist_quilt
        nio_tasks_per_group = 0,
        nio_groups = 2,
        /

