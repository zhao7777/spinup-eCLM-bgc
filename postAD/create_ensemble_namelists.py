#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
from datetime import datetime
from datetime import timedelta
import netCDF4 as nc

def write_modelio(mod, rundir, iens):
    nl_string = ''
    nl_string += '&modelio\n'
    nl_string += '  diri = "UNUSED"\n'
    nl_string += f'  diro = "{rundir}/logs"\n'
    if iens is None:
        nl_string += f'  logfile = "{mod}.log"\n'
    else:
        nl_string += f'  logfile = "{mod}_{str(iens).zfill(4)}.log"\n'
    nl_string += '/\n'
    nl_string += '&pio_inparm\n'
    nl_string += '  pio_netcdf_format = "64bit_offset"\n'
    nl_string += '  pio_numiotasks = -99\n'
    nl_string += '  pio_rearranger = 1\n'
    nl_string += '  pio_root = 1\n'
    nl_string += '  pio_stride = 1\n'
    nl_string += '  pio_typename = "netcdf"\n'
    nl_string += '/'
    if iens is None:
        file_path = f"{rundir}/{mod}_modelio.nml"
    else:
        file_path = f"{rundir}/{mod}_modelio.nml_{str(iens).zfill(4)}"

    with open(file_path, "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)


def write_cpl_modelio(rundir):
    nl_string = ''
    nl_string += '&modelio\n'
    nl_string += '  diri = "UNUSED"\n'
    nl_string += f'  diro = "{rundir}/logs"\n'  
    nl_string += '  logfile = "cpl.log"\n'
    nl_string += '/\n'
    nl_string += '&pio_inparm\n'
    nl_string += '  pio_netcdf_format = "64bit_offset"\n'
    nl_string += '  pio_numiotasks = -99\n'
    nl_string += '  pio_rearranger = 1\n'
    nl_string += '  pio_root = 1\n'
    nl_string += '  pio_stride = 1\n'
    nl_string += '  pio_typename = "netcdf"\n'
    nl_string += '/'
    with open(rundir + "/cpl_modelio.nml", "w+") as nl:
        nl_string = nl_string.replace('"', "'")        
        nl.write(nl_string)

def write_mosart_in(iens, rundir):
    nl_string = ''
    nl_string += '&mosart_inparm\n'
    nl_string += '  bypass_routing_option = "direct_in_place"\n'
    nl_string += '  coupling_period = 10800\n'
    nl_string += '  decomp_option = "roundrobin"\n'
    nl_string += '  delt_mosart = 3600\n'
    nl_string += '  do_rtm = .false.\n'
    nl_string += '  do_rtmflood = .false.\n'
    nl_string += '  finidat_rtm = ""\n'
    nl_string += '  frivinp_rtm = "/p/project1/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm/"\n'
    nl_string += '  ice_runoff = .true.\n'
    nl_string += '  qgwl_runoff_option = "threshold"\n'
    nl_string += '  rtmhist_fexcl1 = ""\n'
    nl_string += '  rtmhist_fexcl2 = ""\n'
    nl_string += '  rtmhist_fexcl3 = ""\n'
    nl_string += '  rtmhist_fincl1 = ""\n'
    nl_string += '  rtmhist_fincl2 = ""\n'
    nl_string += '  rtmhist_fincl3 = ""\n'
    nl_string += '  rtmhist_mfilt = 1\n'
    nl_string += '  rtmhist_ndens = 1\n'
    nl_string += '  rtmhist_nhtfrq = 0\n'
    nl_string += '  smat_option = "Xonly"\n'
    nl_string += '/'
    if iens is None:
        file_path = f"{rundir}/mosart_in"
    else:
        file_path = f"{rundir}/mosart_in_{str(iens).zfill(4)}"

    with open(file_path, "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)

        
def write_datm_in(syear, eyear, iens, rundir):
    nl_string = ''
    nl_string += '&datm_nml\n'
    nl_string += '  decomp = "1d"\n'
    nl_string += '  factorfn = "null"\n'
    nl_string += '  force_prognostic_true = .false.\n'
    nl_string += '  iradsw = 1\n'
    nl_string += '  presaero = .true.\n'
    nl_string += '  restfilm = "undefined"\n'
    nl_string += '  restfils = "undefined"\n'
    nl_string += '  wiso_datm = .false.\n'
    nl_string += '/\n'
    nl_string += '&shr_strdata_nml\n'
    nl_string += '  datamode = "CLMNCEP"\n'
    nl_string += '  domainfile = "/p/project1/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm/domain.lnd.EUR-11_EUR-11.230216_mask.nc"\n'
    if iens is not None:
        nl_string += '  dtlimit = 5.1, 5.1, 5.1, 1.5, 1.5, 5.1, 5.1, 5.1\n'
        nl_string += '  fillalgo = "nn", "nn", "nn", "nn", "nn"\n'
        nl_string += '  fillmask = "nomask", "nomask", "nomask", "nomask", "nomask"\n'
        nl_string += '  fillread = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  fillwrite = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  mapalgo = "nn", "nn", "nn", "bilinear", "bilinear"\n'
        nl_string += '  mapmask = "nomask", "nomask", "nomask", "nomask", "nomask"\n'
        nl_string += '  mapread = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  mapwrite = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  readmode = "single", "single", "single", "single", "single"\n'   
        nl_string += (f'  streams = "datm.streams.ERA5.Solar.stream_' + str(iens).zfill(4) + f' {syear} {syear} {eyear}",\n'
                      f'            "datm.streams.ERA5.Precip.stream_' + str(iens).zfill(4) + f' {syear} {syear} {eyear}",\n'
                      f'            "datm.streams.ERA5.TPQW.stream_' + str(iens).zfill(4) + f' {syear} {syear} {eyear}",\n'
                      '            "datm.streams.presaero.clim_2000 1 2000 2000",\n'
                      '            "datm.streams.topo.observed 1 1 1",\n'
                      f'            "datm.streams.ERA5.Solar_noise.stream_' + str(iens).zfill(4) + f' {syear} {syear} {eyear}",\n'
                      f'            "datm.streams.ERA5.Precip_noise.stream_' + str(iens).zfill(4) + f' {syear} {syear} {eyear}",\n'
                      f'            "datm.streams.ERA5.TPQW_noise.stream_' + str(iens).zfill(4) + f' {syear} {syear} {eyear}",\n')
    else:
        nl_string += '  dtlimit = 5.1, 5.1, 5.1, 1.5, 1.5, \n'
        nl_string += '  fillalgo = "nn", "nn", "nn", "nn", "nn"\n'
        nl_string += '  fillmask = "nomask", "nomask", "nomask", "nomask", "nomask"\n'
        nl_string += '  fillread = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  fillwrite = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  mapalgo = "nn", "nn", "nn", "bilinear", "bilinear"\n'
        nl_string += '  mapmask = "nomask", "nomask", "nomask", "nomask", "nomask"\n'
        nl_string += '  mapread = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  mapwrite = "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET", "NOT_SET"\n'
        nl_string += '  readmode = "single", "single", "single", "single", "single"\n'
        nl_string += (f'  streams = "datm.streams.ERA5.Solar.stream' + f' {syear} {syear} {eyear}",\n'
                      f'            "datm.streams.ERA5.Precip.stream' + f' {syear} {syear} {eyear}",\n'
                      f'            "datm.streams.ERA5.TPQW.stream' + f' {syear} {syear} {eyear}",\n'
                      '            "datm.streams.presaero.clim_2000 1 2000 2000",\n'
                      '            "datm.streams.topo.observed 1 1 1",\n')
    if iens is None:
        nl_string += '  taxmode  = "cycle", "cycle", "cycle", "cycle", "cycle"\n'
        nl_string += '  tintalgo = "coszen", "nearest", "linear", "linear", "lower"\n'
        nl_string += '  vectors = "null"\n'
        nl_string += '/'
        file_path = f"{rundir}/datm_in"
    else:
        nl_string += '  taxmode  = "extend", "extend", "extend", "cycle", "cycle", "extend", "extend", "extend"\n'
        nl_string += '  tintalgo = "nearest", "nearest", "linear", "linear", "lower", "nearest", "nearest", "linear"\n'
        nl_string += '  vectors = "null"\n'
        nl_string += '/'
        file_path = f"{rundir}/datm_in_{str(iens).zfill(4)}"
    with open(file_path, "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)

        
def write_drv_flds_in(rundir):
    nl_string = ""
    nl_string += "&megan_emis_nl\n"
    nl_string += ("  megan_factors_file = '/p/project1/cjibg36/jibg3674/shared_DA/common_eclm_files/megan21_emis_factors_78pft_c20161108.nc'\n")
    nl_string += "  megan_specifier = 'ISOP = isoprene',\n"
    nl_string += "      'C10H16 = pinene_a + carene_3 + thujene_a', 'CH3OH = methanol',\n"
    nl_string += "      'C2H5OH = ethanol', 'CH2O = formaldehyde', 'CH3CHO = acetaldehyde',\n"
    nl_string += "      'CH3COOH = acetic_acid', 'CH3COCH3 = acetone'\n"
    nl_string += "/"

    with open(rundir + "/drv_flds_in", "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)


def write_drv_in(rundir, sdate, edate, ntasks, prefix=None):

    ntasks = str(ntasks)
    nl_string = ""
    nl_string += '&cime_driver_inst\n'
    nl_string += '  ninst_driver = 1\n'
    nl_string += '/\n'
    nl_string += '&cime_pes\n'
    nl_string += '  atm_layout = "concurrent"\n'
    nl_string += f'  atm_ntasks = {ntasks}\n'
    nl_string += '  atm_nthreads = 1\n'
    nl_string += '  atm_pestride = 1\n'
    nl_string += '  atm_rootpe = 0\n'
    nl_string += f'  cpl_ntasks = {ntasks}\n'
    nl_string += '  cpl_nthreads = 1\n'
    nl_string += '  cpl_pestride = 1\n'
    nl_string += '  cpl_rootpe = 0\n'
    nl_string += '  esp_layout = "concurrent"\n'
    nl_string += f'  esp_ntasks = {ntasks}\n'
    nl_string += '  esp_nthreads = 1\n'
    nl_string += '  esp_pestride = 1\n'
    nl_string += '  esp_rootpe = 0\n'
    nl_string += '  glc_layout = "concurrent"\n'
    nl_string += f'  glc_ntasks = {ntasks}\n'
    nl_string += '  glc_nthreads = 1\n'
    nl_string += '  glc_pestride = 1\n'
    nl_string += '  glc_rootpe = 0\n'
    nl_string += '  ice_layout = "concurrent"\n'
    nl_string += f'  ice_ntasks = {ntasks}\n'
    nl_string += '  ice_nthreads = 1\n'
    nl_string += '  ice_pestride = 1\n'
    nl_string += '  ice_rootpe = 0\n'
    nl_string += '  lnd_layout = "concurrent"\n'
    nl_string += f'  lnd_ntasks = {ntasks}\n'
    nl_string += '  lnd_nthreads = 1\n'
    nl_string += '  lnd_pestride = 1\n'
    nl_string += '  lnd_rootpe = 0\n'
    nl_string += '  ocn_layout = "concurrent"\n'
    nl_string += f'  ocn_ntasks = {ntasks}\n'
    nl_string += '  ocn_nthreads = 1\n'
    nl_string += '  ocn_pestride = 1\n'
    nl_string += '  ocn_rootpe = 0\n'
    nl_string += '  rof_layout = "concurrent"\n'
    nl_string += f'  rof_ntasks = {ntasks}\n'
    nl_string += '  rof_nthreads = 1\n'
    nl_string += '  rof_pestride = 1\n'
    nl_string += '  rof_rootpe = 0\n'
    nl_string += '  wav_layout = "concurrent"\n'
    nl_string += f'  wav_ntasks = {ntasks}\n'
    nl_string += '  wav_nthreads = 1\n'
    nl_string += '  wav_pestride = 1\n'
    nl_string += '  wav_rootpe = 0\n'
    nl_string += '/\n'
    nl_string += '&esmf_inparm\n'
    nl_string += '  esmf_logfile_kind = "ESMF_LOGKIND_NONE"\n'
    nl_string += '/\n'
    nl_string += '&papi_inparm\n'
    nl_string += '  papi_ctr1_str = "PAPI_FP_OPS"\n'
    nl_string += '  papi_ctr2_str = "PAPI_NO_CTR"\n'
    nl_string += '  papi_ctr3_str = "PAPI_NO_CTR"\n'
    nl_string += '  papi_ctr4_str = "PAPI_NO_CTR"\n'
    nl_string += '/\n'
    nl_string += '&pio_default_inparm\n'
    nl_string += '  pio_async_interface = .false.\n'
    nl_string += '  pio_blocksize = -1\n'
    nl_string += '  pio_buffer_size_limit = -1\n'
    nl_string += '  pio_debug_level = 0\n'
    nl_string += '  pio_rearr_comm_enable_hs_comp2io = .true.\n'
    nl_string += '  pio_rearr_comm_enable_hs_io2comp = .false.\n'
    nl_string += '  pio_rearr_comm_enable_isend_comp2io = .false.\n'
    nl_string += '  pio_rearr_comm_enable_isend_io2comp = .true.\n'
    nl_string += '  pio_rearr_comm_fcd = "2denable"\n'
    nl_string += '  pio_rearr_comm_max_pend_req_comp2io = 0\n'
    nl_string += '  pio_rearr_comm_max_pend_req_io2comp = 64\n'
    nl_string += '  pio_rearr_comm_type = "p2p"\n'
    nl_string += '/\n'
    nl_string += '&prof_inparm\n'
    nl_string += '  profile_add_detail = .false.\n'
    nl_string += '  profile_barrier = .false.\n'
    nl_string += '  profile_depth_limit = 4\n'
    nl_string += '  profile_detail_limit = 2\n'
    nl_string += '  profile_disable = .false.\n'
    nl_string += '  profile_global_stats = .true.\n'
    nl_string += '  profile_outpe_num = 1\n'
    nl_string += '  profile_outpe_stride = 0\n'
    nl_string += '  profile_ovhd_measurement = .false.\n'
    nl_string += '  profile_papi_enable = .false.\n'
    nl_string += '  profile_single_file = .false.\n'
    nl_string += '  profile_timer = 4\n'
    nl_string += '/\n'
    nl_string += '&seq_cplflds_inparm\n'
    nl_string += '  flds_bgc_oi = .false.\n'
    nl_string += '  flds_co2_dmsa = .false.\n'
    nl_string += '  flds_co2a = .false.\n'
    nl_string += '  flds_co2b = .false.\n'
    nl_string += '  flds_co2c = .false.\n'
    nl_string += '  flds_wiso = .false.\n'
    nl_string += '  glc_nec = 10\n'
    nl_string += '  ice_ncat = 1\n'
    nl_string += '  nan_check_component_fields = .true.\n'
    nl_string += '  seq_flds_i2o_per_cat = .false.\n'
    nl_string += '/\n'
    nl_string += '&seq_cplflds_userspec\n'
    nl_string += '  cplflds_custom = ""\n'
    nl_string += '/\n'
    nl_string += '&seq_flux_mct_inparm\n'
    nl_string += '  seq_flux_atmocn_minwind = 0.5\n'
    nl_string += '  seq_flux_mct_albdif = 0.06\n'
    nl_string += '  seq_flux_mct_albdir = 0.07\n'
    nl_string += '/\n'
    nl_string += '&seq_infodata_inparm\n'
    nl_string += '  aoflux_grid = "ocn"\n'
    nl_string += '  aqua_planet = .false.\n'
    nl_string += '  aqua_planet_sst = 1\n'
    nl_string += '  atm_gnam = "CLM_USRDAT"\n'
    nl_string += '  bfbflag = .false.\n'
    nl_string += '  brnch_retain_casename = .false.\n'
    nl_string += '  budget_ann = 1 \n'
    nl_string += '  budget_daily = 0\n'
    nl_string += '  budget_inst = 0\n'
    nl_string += '  budget_ltann = 1\n'
    nl_string += '  budget_ltend = 0\n'
    nl_string += '  budget_month = 1\n'
    nl_string += '  case_desc = "UNSET"\n'
    if prefix is None:
        nl_string += '  case_name = "clmoas"\n'
    else:
        nl_string += f'  case_name = "{prefix}"\n'
    nl_string += '  cime_model = "cesm"\n'
    nl_string += '  coldair_outbreak_mod = .true.\n'
    nl_string += '  cpl_decomp = 0\n'
    nl_string += '  cpl_seq_option = "CESM1_MOD"\n'
    nl_string += '  do_budgets = .false.\n'
    nl_string += '  do_histinit = .false.\n'
    nl_string += '  drv_threading = .false.\n'
    nl_string += '  eps_aarea = 9e-07\n'
    nl_string += '  eps_agrid = 1e-12\n'
    nl_string += '  eps_amask = 1e-13\n'
    nl_string += '  eps_frac = 0.01\n'
    nl_string += '  eps_oarea = 0.1\n'
    nl_string += '  eps_ogrid = 0.01\n'
    nl_string += '  eps_omask = 1e-06\n'
    nl_string += '  flux_albav = .false.\n'
    nl_string += '  flux_convergence = 0.01\n'
    nl_string += '  flux_diurnal = .false.\n'
    nl_string += '  flux_epbal = "off"\n'
    nl_string += '  flux_max_iteration = 5\n'
    nl_string += '  force_stop_at = "month"\n'
    nl_string += '  glc_gnam = "null"\n'
    nl_string += '  glc_renormalize_smb = "on_if_glc_coupled_fluxes"\n'
    nl_string += '  gust_fac = 0.0\n'
    nl_string += '  histaux_a2x = .false.\n'
    nl_string += '  histaux_a2x1hr = .false.\n'
    nl_string += '  histaux_a2x1hri = .false.\n'
    nl_string += '  histaux_a2x24hr = .false.\n'
    nl_string += '  histaux_a2x3hr = .false.\n'
    nl_string += '  histaux_a2x3hrp = .false.\n'
    nl_string += '  histaux_double_precision = .false.\n'
    nl_string += '  histaux_l2x = .false.\n'
    nl_string += '  histaux_l2x1yrg = .false.\n'
    nl_string += '  histaux_r2x = .false.\n'
    nl_string += '  histavg_atm = .true.\n'
    nl_string += '  histavg_glc = .true.\n'
    nl_string += '  histavg_ice = .true.\n'
    nl_string += '  histavg_lnd = .true.\n'
    nl_string += '  histavg_ocn = .true.\n'
    nl_string += '  histavg_rof = .true.\n'
    nl_string += '  histavg_wav = .true.\n'
    nl_string += '  histavg_xao = .true.\n'
    nl_string += '  hostname = "juwels"\n'
    nl_string += '  ice_gnam = "null"\n'
    nl_string += '  info_debug = 1\n'
    nl_string += '  lnd_gnam = "CLM_USRDAT"\n'
    nl_string += '  logfilepostfix = ".log"\n'
    nl_string += '  max_cplstep_time = 0.0\n'
    nl_string += '  mct_usealltoall = .false.\n'
    nl_string += '  mct_usevector = .false.\n'
    nl_string += '  model_doi_url = "https://doi.org/10.5065/D67H1H0V"\n'
    nl_string += '  model_version = "release-clm5.0.34-2-ga2989b0"\n'
    nl_string += '  ocn_gnam = "null"\n'
    nl_string += '  orb_eccen = 1e+36\n'
    nl_string += '  orb_iyear = 2000\n'
    nl_string += '  orb_iyear_align = 2000\n'
    nl_string += '  orb_mode = "fixed_year"\n'
    nl_string += '  orb_mvelp = 1e+36\n'
    nl_string += '  orb_obliq = 1e+36\n'
    nl_string += f'  outpathroot = "{rundir}/out/"\n'
    nl_string += '  reprosum_diffmax = -1.0e-8\n'
    nl_string += '  reprosum_recompute = .false.\n'
    nl_string += '  reprosum_use_ddpdd = .false.\n'
    nl_string += '  restart_file = "str_undefined"\n'
    nl_string += '  rof_gnam = "null"\n'
    nl_string += '  run_barriers = .false.\n'
    nl_string += '  scmlat = -999.\n'
    nl_string += '  scmlon = -999.\n'
    nl_string += '  shr_map_dopole = .true.\n'
    nl_string += '  single_column = .false.\n'
    nl_string += '  start_type = "startup"\n'
    nl_string += f'  tchkpt_dir = "{rundir}/timing/checkpoints"\n'
    nl_string += '  tfreeze_option = "mushy"\n'
    nl_string += f'  timing_dir = "{rundir}/timing"\n'
    nl_string += '  username = "zhao2"\n'
    nl_string += '  vect_map = "cart3d"\n'
    nl_string += '  wall_time_limit = -1.0\n'
    nl_string += '  wav_gnam = "null"\n'
    nl_string += '  wv_sat_scheme = "GoffGratch"\n'
    nl_string += '  wv_sat_table_spacing = 1.0\n'
    nl_string += '  wv_sat_transition_start = 20.0\n'
    nl_string += '  wv_sat_use_tables = .false.\n'
    nl_string += '/\n'
    nl_string += '&seq_timemgr_inparm\n'
    nl_string += '  atm_cpl_dt = 3600\n'
    nl_string += '  atm_cpl_offset = 0\n'
    nl_string += '  barrier_n = 1\n'
    nl_string += '  barrier_option = "ndays"\n'
    nl_string += '  barrier_ymd = -999\n'
    nl_string += '  calendar = "NO_LEAP"\n'
    nl_string += '  data_assimilation_atm = .false.\n'
    nl_string += '  data_assimilation_cpl = .false.\n'
    nl_string += '  data_assimilation_glc = .false.\n'
    nl_string += '  data_assimilation_ice = .false.\n'
    nl_string += '  data_assimilation_lnd = .false.\n'
    nl_string += '  data_assimilation_ocn = .false.\n'
    nl_string += '  data_assimilation_rof = .false.\n'
    nl_string += '  data_assimilation_wav = .false.\n'
    nl_string += '  end_restart = .false.\n'
    nl_string += '  esp_cpl_offset = 0\n'
    nl_string += '  esp_run_on_pause = .true.\n'
    nl_string += '  glc_avg_period = "yearly"\n'
    nl_string += '  glc_cpl_dt = 3600\n'
    nl_string += '  glc_cpl_offset = 0\n'
    nl_string += '  histavg_n = -999\n'
    nl_string += '  histavg_option = "never"\n'
    nl_string += '  histavg_ymd = -999\n'
    nl_string += '  history_n = -999\n'
    nl_string += '  history_option = "never"\n'
    nl_string += '  history_ymd = -999\n'
    nl_string += '  ice_cpl_dt = 3600\n'
    nl_string += '  ice_cpl_offset = 0\n'
    nl_string += '  lnd_cpl_dt = 3600\n'
    nl_string += '  lnd_cpl_offset = 0\n'
    nl_string += '  ocn_cpl_dt = 3600\n'
    nl_string += '  ocn_cpl_offset = 0\n'
    nl_string += '  pause_active_atm = .false.\n'
    nl_string += '  pause_active_cpl = .false.\n'
    nl_string += '  pause_active_glc = .false.\n'
    nl_string += '  pause_active_ice = .false.\n'
    nl_string += '  pause_active_lnd = .false.\n'
    nl_string += '  pause_active_ocn = .false.\n'
    nl_string += '  pause_active_rof = .false.\n'
    nl_string += '  pause_active_wav = .false.\n'
    nl_string += '  pause_n = 0\n'
    nl_string += '  pause_option = "never"\n'
    # nl_string += '  restart_n = -1\n'
    # nl_string += '  restart_option = "date"\n'
    # nl_string += f'  restart_ymd = {sdate}\n'
    nl_string += '  restart_n = 1\n'
    nl_string += '  restart_option = "nyears"\n'
    nl_string += '  restart_ymd = -1\n'
    nl_string += '  rof_cpl_dt = 10800\n'
    nl_string += '  start_tod = 0\n'
    nl_string += f'  start_ymd = {sdate}\n'
    # nl_string += '  stop_n = -1\n'
    # nl_string += '  stop_option = "date"\n'
    nl_string += '  stop_n = 120\n'
    nl_string += '  stop_option = "nmonths"\n'
    # nl_string += f'  stop_ymd = {edate}\n'
    nl_string += '  stop_ymd = -1\n'
    nl_string += '  tprof_n = -999\n'
    nl_string += '  tprof_option = "never"\n'
    nl_string += '  tprof_ymd = -999\n'
    nl_string += '  wav_cpl_dt = 3600\n'
    nl_string += '  wav_cpl_offset = 0\n'
    nl_string += '/\n'

    with open(rundir + "/drv_in", "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)
        
        
def write_lnd_in(rundir, iens, restart_from=None):
    nl_string = "" 
    nl_string += "&clm_inparm\n"
    nl_string += " albice = 0.5,0.3\n"
    # nl_string += " co2_ppmv = 367.0\n"
    nl_string += " co2_ppmv = 315.98\n"
    nl_string += " co2_type = 'constant'\n"
    nl_string += " create_crop_landunit = .true.\n"
    nl_string += " dtime = 3600\n"
    nl_string += " fatmlndfrc = '/p/project1/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm/domain.lnd.EUR-11_EUR-11.230216_mask.nc'\n"
    # nl_string += f" finidat = '{finidat}'\n"
    if iens is not None:
        nl_string += " finidat = '/p/scratch/detectrea/zhao2/CLM5_DA/DA_eCLM_cordex_444x432_1y_iter5/20190101-20191230/i005/R" + str(iens + 1).zfill(3) + "/run_20181003-20190101/EU11.clm2.r.2019-01-01-00000.nc'\n"
    else:
        # iens_rn = 0
        # nl_string += " finidat = '/p/scratch/detectrea/zhao2/CLM5_DA/DA_eCLM_cordex_444x432_1y_iter5/20190101-20191230/i005/R" + str(iens_rn + 1).zfill(3) + "/run_20181003-20190101/EU11.clm2.r.2019-01-01-00000.nc'\n"
        nl_string += f" finidat = '{restart_from}'\n"
    nl_string += " fsnowaging = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/snicar_drdt_bst_fit_60_c070416.nc'\n"
    nl_string += " fsnowoptics = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/snicar_optics_5bnd_c090915.nc'\n"
    if iens is not None:
        nl_string += " fsurdat = '/p/scratch/detectrea/zhao2/CLM5_DA/DA_eCLM_cordex_444x432_1y_iter5/20190101-20191230/i005/R" + str(iens + 1).zfill(3) + "/surfdata_EUR-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230808_GLC2000.nc'\n"
    else:
        iens_rn = 0
        # nl_string += " fsurdat = '/p/scratch/detectrea/zhao2/CLM5_DA/DA_eCLM_cordex_444x432_1y_iter5/input_clm/surfdata_EUR-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230808_GLC2000.nc'\n"
        nl_string += " fsurdat = '/p/project/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm/surfdata_EUR-11_hist_16pfts_Irrig_CMIP6_simyr2000_c230808_GLC2000.nc'\n"
    nl_string += " glc_do_dynglacier = .false.\n"
    nl_string += " glc_snow_persistence_max_days = 0\n"
    nl_string += " h2osno_max = 10000.0\n"
    nl_string += " hist_empty_htapes = .true.\n"
    nl_string += " hist_dov2xy = .true.\n"
    # nl_string += " hist_fincl1 = 'QFLX_EVAP_TOT','SOILLIQ', 'SOILICE', 'H2OSOI', 'ERRH2O', 'QDRAI', 'QRUNOFF', 'QOVER','TWS', 'FSNO', 'H2OSNO', 'RAIN'\n"
    # nl_string += " hist_mfilt = 365\n"
    # nl_string += " hist_nhtfrq = -24\n"
    nl_string += " hist_fincl1 = 'GPP','TLAI','TOTECOSYSC','TOTECOSYSN','TOTSOMC','TOTVEGC','TWS'\n"
    nl_string += " hist_mfilt = 12\n"
    nl_string += " hist_nhtfrq = 0\n"    
    nl_string += " int_snow_max = 2000\n"
    nl_string += " irrigate = .true.\n"
    nl_string += " maxpatch_glcmec = 10.0\n"
    nl_string += " maxpatch_pft = 79\n"
    nl_string += " n_melt_glcmec = 10.0\n"
    nl_string += " nlevsno = 12\n"
    nl_string += " nsegspc = 35\n"
    if iens is not None:
        nl_string += " paramfile = '/p/scratch/detectrea/zhao2/CLM5_DA/DA_eCLM_cordex_444x432_1y_iter5/20190101-20191230/i005/R" + str(iens + 1).zfill(3) + "/clm5_params.c171117.nc'\n"
    else:
        # nl_string += " paramfile = /p/scratch/detectrea/zhao2/CLM5_DA/DA_eCLM_cordex_444x432_1y_iter5/input_clm/clm5_params.c171117.nc'\n"
        nl_string += " paramfile = '/p/project/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm/clm5_params.c171117.nc'\n"
    nl_string += " run_zero_weight_urban = .false.\n"
    nl_string += " soil_layerstruct = '20SL_8.5m'\n"
    nl_string += " spinup_state = 0\n"
    nl_string += " suplnitro = 'NONE'\n"
    nl_string += " use_bedrock = .true.\n"
    nl_string += " use_century_decomp = .true.\n"
    nl_string += " use_cn = .true.\n"
    nl_string += " use_crop = .false.\n"
    nl_string += " use_dynroot = .false.\n"
    nl_string += " use_fates = .false.\n"
    nl_string += " use_fertilizer = .true.\n"
    nl_string += " use_flexiblecn = .true.\n"
    nl_string += " use_fun = .true.\n"
    nl_string += " use_grainproduct = .true.\n"
    nl_string += " use_hydrstress = .true.\n"
    # nl_string += " use_init_interp = .false.\n"
    nl_string += " use_lai_streams = .false.\n"
    nl_string += " use_lch4 = .true.\n"
    nl_string += " use_luna = .true.\n"
    nl_string += " use_nguardrail = .true.\n"
    nl_string += " use_nitrif_denitrif = .true.\n"
    nl_string += " use_soil_moisture_streams = .false.\n"
    nl_string += " use_vertsoilc = .true.\n"
    nl_string += "/\n"
    nl_string += "&ndepdyn_nml\n"
    nl_string += " ndep_taxmode = 'cycle'\n" 
    nl_string += " ndep_varlist = 'NDEP_month'\n" 
    nl_string += " ndepmapalgo = 'bilinear'\n" 
    nl_string += " stream_fldfilename_ndep = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/fndep_clm_hist_b.e21.BWHIST.f09_g17.CMIP6-historical-WACCM.ensmean_1849-2015_monthly_0.9x1.25_c180926.nc'\n"
    nl_string += " stream_year_first_ndep = 2000\n"
    nl_string += " stream_year_last_ndep = 2000\n"
    nl_string += "/\n"
    nl_string += "&popd_streams\n"
    nl_string += " popdensmapalgo = 'bilinear'\n" 
    nl_string += " stream_fldfilename_popdens = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/clmforc.Li_2017_HYDEv3.2_CMIP6_hdm_0.5x0.5_AVHRR_simyr1850-2016_c180202.nc'\n"
    nl_string += " stream_year_first_popdens = 2000\n"
    nl_string += " stream_year_last_popdens = 2000\n"    
    nl_string += "/\n"
    nl_string += "&urbantv_streams\n"
    nl_string += " stream_fldfilename_urbantv = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/CLM50_tbuildmax_Oleson_2016_0.9x1.25_simyr1849-2106_c160923.nc'\n"
    nl_string += " stream_year_first_urbantv = 2000\n"
    nl_string += " stream_year_last_urbantv = 2000\n"
    nl_string += " urbantvmapalgo = 'nn'\n"
    nl_string += "/\n"
    nl_string += "&light_streams\n"
    nl_string += " lightngmapalgo = 'bilinear'\n" 
    nl_string += " stream_fldfilename_lightng = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/clmforc.Li_2012_climo1995-2011.T62.lnfm_Total_c140423.nc'\n"
    nl_string += " stream_year_first_lightng = 0001\n"
    nl_string += " stream_year_last_lightng = 0001\n"      
    nl_string += "/\n"
    nl_string += "&soil_moisture_streams\n"
    nl_string += "/\n"
    nl_string += "&lai_streams\n"
    nl_string += "/\n"
    nl_string += "&atm2lnd_inparm\n"
    nl_string += " glcmec_downscale_longwave = .true.\n"
    nl_string += " lapse_rate = 0.006\n"
    nl_string += " lapse_rate_longwave = 0.032\n"
    nl_string += " longwave_downscaling_limit = 0.5\n"
    nl_string += " precip_repartition_glc_all_rain_t = 0.\n"
    nl_string += " precip_repartition_glc_all_snow_t = -2.\n"
    nl_string += " precip_repartition_nonglc_all_rain_t = 2.\n"
    nl_string += " precip_repartition_nonglc_all_snow_t = 0.\n"
    nl_string += " repartition_rain_snow = .true.\n"
    nl_string += "/\n"
    nl_string += "&lnd2atm_inparm\n"
    nl_string += " melt_non_icesheet_ice_runoff = .true.\n"
    nl_string += "/\n"
    nl_string += "&clm_canopyhydrology_inparm\n"
    nl_string += " interception_fraction = 1.0\n"
    nl_string += " maximum_leaf_wetted_fraction = 0.05\n"
    nl_string += " snowveg_flag = 'ON_RAD'\n"
    nl_string += " use_clm5_fpi = .true.\n"
    nl_string += "/\n"
    nl_string += "&cnphenology\n"
    nl_string += " initial_seed_at_planting = 3.d00\n"    
    nl_string += "/\n"
    nl_string += "&clm_soilhydrology_inparm\n"
    nl_string += "/\n"
    nl_string += "&dynamic_subgrid\n"
    nl_string += "/\n"
    nl_string += "&cnvegcarbonstate\n"
    nl_string += " initial_vegc = 100.d00\n"     
    nl_string += "/\n"
    nl_string += "&finidat_consistency_checks\n"
    nl_string += "/\n"
    nl_string += "&dynpft_consistency_checks\n"
    nl_string += "/\n"
    nl_string += "&clm_initinterp_inparm\n"
    nl_string += " init_interp_method = 'general'\n"
    nl_string += "/\n"
    nl_string += "&century_soilbgcdecompcascade\n"
    nl_string += " initial_cstocks = 200.0d00, 200.0d00, 200.0d00\n"     
    nl_string += " initial_cstocks_depth = 1.50d00\n"         
    nl_string += "/\n"
    nl_string += "&soilhydrology_inparm\n"
    nl_string += " baseflow_scalar = 0.001d00\n"
    nl_string += "/\n"
    nl_string += "&luna\n"
    nl_string += " jmaxb1 = 0.093563\n"
    nl_string += "/\n"
    nl_string += "&friction_velocity\n"
    nl_string += " zetamaxstable = 0.5\n"
    nl_string += "/\n"
    nl_string += "&mineral_nitrogen_dynamics\n"
    nl_string += "/\n"
    nl_string += "&soilwater_movement_inparm\n"
    nl_string += " dtmin = 60.0\n"
    nl_string += " expensive = 42\n"
    nl_string += " flux_calculation = 1\n"
    nl_string += " inexpensive = 1\n"
    nl_string += " lower_boundary_condition = 2\n"
    nl_string += " soilwater_movement_method = 1\n"
    nl_string += " upper_boundary_condition = 1\n"
    nl_string += " verysmall = 1.e-8\n"
    nl_string += " xtolerlower = 1.e-2\n"
    nl_string += " xtolerupper = 1.e-1\n"
    nl_string += "/\n"
    nl_string += "&rooting_profile_inparm\n"
    nl_string += " rooting_profile_method_carbon = 1\n"
    nl_string += " rooting_profile_method_water = 1\n"
    nl_string += "/\n"
    nl_string += "&soil_resis_inparm\n"
    nl_string += " soil_resis_method = 1\n"
    nl_string += "/\n"
    nl_string += "&bgc_shared\n"
    nl_string += " constrain_stress_deciduous_onset = .true.\n"
    nl_string += " decomp_depth_efolding = 10.0\n"    
    nl_string += "/\n"
    nl_string += "&canopyfluxes_inparm\n"
    nl_string += " use_undercanopy_stability = .false.\n"
    nl_string += "/\n"
    nl_string += "&aerosol\n"
    nl_string += " fresh_snw_rds_max = 204.526\n"
    nl_string += "/\n"
    nl_string += "&clmu_inparm\n"
    nl_string += " building_temp_method = 1\n"
    nl_string += " urban_hac = 'ON_WASTEHEAT'\n"
    nl_string += " urban_traffic = .false.\n"
    nl_string += "/\n"
    nl_string += "&clm_soilstate_inparm\n"
    nl_string += " organic_frac_squared = .false.\n"
    nl_string += "/\n"
    nl_string += "&clm_nitrogen\n"
    nl_string += "  carbon_resp_opt = 0\n"
    nl_string += "  cn_evergreen_phenology_opt = 1\n"
    nl_string += "  cn_partition_opt  = 1\n"
    nl_string += "  cn_residual_opt = 1\n"
    nl_string += "  cnratio_floating  = .true.\n"
    nl_string += "  downreg_opt = .false.\n"
    nl_string += "  lnc_opt = .true.\n"
    nl_string += "  mm_nuptake_opt = .true.\n"
    nl_string += "  nscalar_opt = .true.\n"
    nl_string += "  plant_ndemand_opt = 3\n"
    nl_string += "  reduce_dayl_factor = .false.\n"
    nl_string += "  substrate_term_opt = .true.\n"
    nl_string += "  temp_scalar_opt = .true.\n"
    nl_string += "  vcmax_opt = 3\n"
    nl_string += " lnc_opt = .false.\n"
    nl_string += "/\n"
    nl_string += "&clm_snowhydrology_inparm\n"
    nl_string += " lotmp_snowdensity_method = 'Slater2017'\n"
    nl_string += " reset_snow = .false.\n"
    nl_string += " reset_snow_glc = .false.\n"
    nl_string += " reset_snow_glc_ela = 1.e9\n"
    nl_string += " snow_overburden_compaction_method = 'Vionnet2012'\n"
    nl_string += " upplim_destruct_metamorph = 175.0\n"
    nl_string += " wind_dependent_snow_density = .true.\n"
    nl_string += "/\n"
    nl_string += "&cnprecision_inparm\n"
    nl_string += "  cnegcrit = -6.d+1\n"
    nl_string += "  ncrit = 1.d-9\n"
    nl_string += "  nnegcrit = -6.d+0\n"
    nl_string += "/\n"
    nl_string += "&clm_glacier_behavior\n"
    nl_string += " glacier_region_behavior = 'single_at_atm_topo','virtual','virtual','multiple'\n"
    nl_string += " glacier_region_ice_runoff_behavior = 'melted','melted','remains_ice','remains_ice'\n"
    nl_string += " glacier_region_melt_behavior = 'remains_in_place','replaced_by_ice','replaced_by_ice','replaced_by_ice'\n"
    nl_string += " glacier_region_rain_to_snow_behavior = 'converted_to_snow','converted_to_snow','converted_to_snow','converted_to_snow'\n"
    nl_string += "/\n"
    nl_string += "&crop\n"
    nl_string += "  baset_latvary_intercept = 12.0d00\n"
    nl_string += "  baset_latvary_slope = 0.4d00\n"
    nl_string += "  baset_mapping = 'varytropicsbylat'\n"
    nl_string += "/\n"
    nl_string += "&irrigation_inparm\n"
    nl_string += " irrig_depth = 0.6\n"
    nl_string += " irrig_length = 14400\n"
    nl_string += " irrig_min_lai = 0.0\n"
    nl_string += " irrig_start_time = 21600\n"
    nl_string += " irrig_target_smp = -3400.\n"
    nl_string += " irrig_threshold_fraction = 1.0\n"
    nl_string += " limit_irrigation_if_rof_enabled = .false.\n"
    nl_string += "/\n"
    nl_string += "&ch4par_in\n"
    nl_string += "  finundation_method = 'TWS_inversion'\n"
    nl_string += "  use_aereoxid_prog = .true.\n"
    nl_string += "/\n"
    nl_string += "&clm_humanindex_inparm\n"
    nl_string += " calc_human_stress_indices = 'FAST'\n"
    nl_string += "/\n"
    nl_string += "&cnmresp_inparm\n"
    nl_string += " br_root = 0.83d-06\n"
    nl_string += "/\n"
    nl_string += "&photosyns_inparm\n"
    nl_string += " leafresp_method = 0\n"
    nl_string += " light_inhibit = .true.\n"
    nl_string += " modifyphoto_and_lmr_forcrop = .true.\n"
    nl_string += " rootstem_acc = .false.\n"
    nl_string += " stomatalcond_method = 'Medlyn2011'\n"
    nl_string += "/\n"
    nl_string += "&cnfire_inparm\n"
    nl_string += "  fire_method = 'li2016crufrc'\n"
    nl_string += "/\n"
    nl_string += "&cn_general\n"
    nl_string += "  dribble_crophrv_xsmrpool_2atm = .false.\n"
    nl_string += "/\n"
    nl_string += "&nitrif_inparm\n"
    nl_string += "/\n"
    nl_string += "&lifire_inparm\n"
    nl_string += "  boreal_peatfire_c = 0.09d-4\n"
    nl_string += "  bt_max = 0.98d00\n"
    nl_string += "  bt_min = 0.85d00\n"
    nl_string += "  cli_scale = 0.033d00\n"
    nl_string += "  cmb_cmplt_fact = 0.5d00, 0.28d00\n"
    nl_string += "  cropfire_a1 = 1.6d-4\n"
    nl_string += "  lfuel = 105.d00\n"
    nl_string += "  non_boreal_peatfire_c = 0.17d-3\n"
    nl_string += "  occur_hi_gdp_tree = 0.33d00\n"
    nl_string += "  pot_hmn_ign_counts_alpha = 0.010d00\n"
    nl_string += "  rh_hgh = 80.0d00\n"
    nl_string += "  rh_low = 30.0d00\n"
    nl_string += "  ufuel = 1050.d00\n"    
    nl_string += "/\n"
    nl_string += "&ch4finundated\n"
    nl_string += "  stream_fldfilename_ch4finundated = '/p/project/cjibg36/jibg3674/shared_DA/common_eclm_files/finundated_inversiondata_0.9x1.25_c170706.nc'\n"
    nl_string += "/\n"
    nl_string += "&clm_canopy_inparm\n"
    nl_string += " leaf_mr_vcm = 0.015\n"
    nl_string += "/\n"
    nl_string += ""
    if iens is None:
        file_path = f"{rundir}/lnd_in"
    else:
        file_path = f"{rundir}/lnd_in_{str(iens).zfill(4)}"
    with open(file_path, "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)


def write_stream_files(rundir, iens, nreal, sdate, edate, stream):
    nl_string = ''
    nl_string += '<?xml version="1.0"?>\n'
    nl_string += '<file id="stream" version="1.0">\n'
    nl_string += '  <dataSource>\n'
    nl_string += '    GENERIC\n'
    nl_string += '  </dataSource>\n'
    nl_string += '  <domainInfo>\n'
    nl_string += '    <variableNames>\n'
    nl_string += '     time    time\n'
    nl_string += '     xc      lon\n'
    nl_string += '     yc      lat\n'
    nl_string += '     area    area\n'
    nl_string += '     mask    mask\n'
    nl_string += '    </variableNames>\n'
    nl_string += '    <filePath>\n'
    nl_string += '     /p/project1/cjibg36/jibg3674/shared_DA/setup_eclm_cordex_444x432_v9/input_clm\n'
    nl_string += '    </filePath>\n'
    nl_string += '    <fileNames>\n'
    nl_string += '     domain.lnd.EUR-11_EUR-11.230216_mask.nc\n'
    nl_string += '    </fileNames>\n'
    nl_string += '  </domainInfo>\n'
    nl_string += '  <fieldInfo>\n'
    nl_string += '    <variableNames>\n'
    if stream == "TPQW":
        # nl_string += '      Tdew     tdew\n'        
        nl_string += '      PSRF     pbot\n'
        nl_string += '      QBOT     shum\n'
        nl_string += '      TBOT     tbot\n'
        nl_string += '      WIND     wind\n'
        nl_string += '      FLDS     lwdn\n'
    elif stream == "Precip":
        nl_string += '      PRECTmms precn\n'
    elif stream == "Solar":
        nl_string += '      FSDS swdn\n'
    elif stream == "TPQW_noise":
        nl_string += '      TBOT     tbot_noise\n'
        nl_string += '      FLDS     lwdn_noise\n'  
    elif stream == "Precip_noise":
        nl_string += '      PRECTmms     precn_noise\n'
    elif stream == "Solar_noise":
        nl_string += '      FSDS     swdn_noise\n'
    nl_string += '    </variableNames>\n'
    nl_string += '    <filePath>\n'
    if stream in ["TPQW", "Precip", "Solar"]:
        nl_string += '      /p/scratch/cjibg36/jibg3674/ERA5_forcing/data_out/\n'
    if stream in ["TPQW_noise", "Precip_noise", "Solar_noise"]:
        nl_string += '      /p/project1/detectc01/clm_inputfiles/forcings/noise/\n'
    nl_string += '    </filePath>\n'
    nl_string += '    <fileNames>\n'
    # sdate = str(sdate)
    # edate = str(edate)
    start_date = datetime.strptime(sdate, "%Y%m%d")
    end_date = datetime.strptime(edate, "%Y%m%d")
    current_date = start_date
    while (current_date.year, current_date.month) <= (end_date.year, end_date.month):
        year_month = current_date.strftime("%Y-%m")  
        nl_string += f'      {year_month}.nc\n'
        if current_date.month == 12:
            current_date = current_date.replace(year=current_date.year + 1, month=1, day=1)
        else:
            current_date = current_date.replace(month=current_date.month + 1, day=1)
    nl_string += '    </fileNames>\n'
    nl_string += '  <offset>\n'
    nl_string += '    0 \n'
    nl_string += '  </offset>\n'
    if stream in ["TPQW_noise", "Precip_noise", "Solar_noise"]:
        nl_string += '  <timeInformation>\n'
        nl_string += '    3\n'
        nl_string += '  </timeInformation>\n'
        nl_string += '  <caseId>\n'
        if iens is None:
            iens_rn = 0
            nl_string += f'    {iens_rn}\n'
        else:
            nl_string += f'    {iens}\n'
            nl_string += '  </caseId>\n'
            nl_string += '  <numEns>\n'
            nl_string += f'    {nreal}\n'
            nl_string += '  </numEns>\n'
    nl_string += '  </fieldInfo>\n'
    nl_string +='</file>\n'
    nl_string +=''
    if iens is None:
        file_path = f"{rundir}/datm.streams.ERA5.{stream}.stream"
    else:
        file_path = f"{rundir}/datm.streams.ERA5.{stream}.stream_" + str(iens).zfill(4)
    with open(file_path, "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)

        
def write_presaero_stream_files(rundir):
    nl_string = ''
    nl_string += '<?xml version="1.0"?>\n'
    nl_string += '<file id="stream" version="1.0">\n'
    nl_string += '  <dataSource>\n'
    nl_string += '    GENERIC\n'
    nl_string += '  </dataSource>\n'
    nl_string += '  <domainInfo>\n'
    nl_string += '    <variableNames>\n'
    nl_string += '      time    time\n'
    nl_string += '      lon     lon\n'
    nl_string += '      lat     lat\n'
    nl_string += '      area    area\n'
    nl_string += '      mask    mask\n'
    nl_string += '    </variableNames>\n'
    nl_string += '    <filePath>\n'
    nl_string += '      /p/project1/cjibg36/jibg3674/shared_DA/common_eclm_files\n'
    nl_string += '    </filePath>\n'
    nl_string += '    <fileNames>\n'
    nl_string += '      aerosoldep_WACCM.ensmean_monthly_hist_1849-2015_0.9x1.25_CMIP6_c180926.nc\n'
    nl_string += '    </fileNames>\n'
    nl_string += '  </domainInfo>\n'
    nl_string += '  <fieldInfo>\n'
    nl_string += '    <variableNames>\n'
    nl_string += '      BCDEPWET   bcphiwet\n'
    nl_string += '      BCPHODRY   bcphodry\n'
    nl_string += '      BCPHIDRY   bcphidry\n'
    nl_string += '      OCDEPWET   ocphiwet\n'
    nl_string += '      OCPHIDRY   ocphidry\n'
    nl_string += '      OCPHODRY   ocphodry\n'
    nl_string += '      DSTX01WD   dstwet1\n'
    nl_string += '      DSTX01DD   dstdry1\n'
    nl_string += '      DSTX02WD   dstwet2\n'
    nl_string += '      DSTX02DD   dstdry2\n'
    nl_string += '      DSTX03WD   dstwet3\n'
    nl_string += '      DSTX03DD   dstdry3\n'
    nl_string += '      DSTX04WD   dstwet4\n'
    nl_string += '      DSTX04DD   dstdry4\n'
    nl_string += '    </variableNames>\n'
    nl_string += '    <filePath>\n'
    nl_string += '      /p/project1/cjibg36/jibg3674/shared_DA/common_eclm_files\n'
    nl_string += '    </filePath>\n'
    nl_string += '    <fileNames>\n'
    nl_string += '      aerosoldep_WACCM.ensmean_monthly_hist_1849-2015_0.9x1.25_CMIP6_c180926.nc\n'
    nl_string += '    </fileNames>\n'
    nl_string += '    <offset>\n'
    nl_string += '      0\n'
    nl_string += '    </offset>\n'
    nl_string += '  </fieldInfo>\n'
    nl_string += '</file>\n'
    
    with open(rundir + "/datm.streams.presaero.clim_2000",  "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)

def write_topo_stream_files(rundir):
    nl_string = ''
    nl_string += '<?xml version="1.0"?>\n'
    nl_string += '<file id="stream" version="1.0">\n'
    nl_string += '  <dataSource>\n'
    nl_string += '    GENERIC\n'
    nl_string += '  </dataSource>\n'
    nl_string += '  <domainInfo>\n'
    nl_string += '    <variableNames>\n'
    nl_string += '      time    time\n'
    nl_string += '      LONGXY  lon\n'
    nl_string += '      LATIXY  lat\n'
    nl_string += '      area    area\n'
    nl_string += '      mask    mask\n'
    nl_string += '    </variableNames>\n'
    nl_string += '    <filePath>\n'
    nl_string += '      /p/project1/cjibg36/jibg3674/shared_DA/common_eclm_files\n'
    nl_string += '    </filePath>\n'
    nl_string += '    <fileNames>\n'
    nl_string += '      topodata_0.9x1.25_USGS_070110_stream_c151201.nc\n'
    nl_string += '    </fileNames>\n'
    nl_string += '  </domainInfo>\n'
    nl_string += '  <fieldInfo>\n'
    nl_string += '    <variableNames>\n'
    nl_string += '      TOPO    topo\n'
    nl_string += '    </variableNames>\n'
    nl_string += '    <filePath>\n'
    nl_string += '      /p/project1/cjibg36/jibg3674/shared_DA/common_eclm_files\n'
    nl_string += '    </filePath>\n'
    nl_string += '    <fileNames>\n'
    nl_string += '      topodata_0.9x1.25_USGS_070110_stream_c151201.nc\n'
    nl_string += '    </fileNames>\n'
    nl_string += '  <offset>\n'
    nl_string += '    0\n'
    nl_string += '  </offset>\n'
    nl_string += '  </fieldInfo>\n'
    nl_string += '</file>\n'
    
    with open(rundir + "/datm.streams.topo.observed",  "w+") as nl:
        nl_string = nl_string.replace('"', "'")
        nl.write(nl_string)

