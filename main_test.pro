pro main_test
;;; test routines
computer = 'I:'
@folders

del_data, '*'
;trange = ['2011-06-21/13:09:23', '2011-06-21/13:13:23']
trange = ['2009-04-07/07:03:23', '2009-04-07/07:08:23']
probe = 'd'
;vdht = vdht_calc(trange, probe_arr = probe, fgs_folder = fgs_folder, efs_folder = efs_folder, maxgap = 10)

thm_load_state, probe = probe, trange = trange, /get_supp
thm_load_fit, probe = probe, trange = trange, coord='dsl', level = 2, data='efs'
thm_load_efi, probe = probe, trange = trange, coord='dsl', level = 2, data='eff_dot0';, /get_supp
thm_load_efi, probe = probe, trange = trange, coord='dsl', level = 1, data='eff_0'
get_data, 'thd_eff_dot0_dsl', t, data
store_data, 'thd_eff_dot0_dsl_xy', data = {x:t, y:data(*, 0:1)}
options, 'thd_eff_dot0_dsl_xy', ysubtitle = '[mV/m]', colors = [2,4]

tplot, ['thd_efs_dsl', 'thd_efs_dot0_dsl', 'thd_eff_dot0_dsl_xy', 'thd_eff_0']

stop
end
