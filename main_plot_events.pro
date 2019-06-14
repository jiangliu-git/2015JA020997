pro main_plot_events
;;; check events
;; plot an example of multi-point observation
thm_init
computer = 'I:'
;computer = '/home/jliu'
@folders

;; some constants for loading data
minutes_load = 1.

method = 'binxbout'
v_use = 'efs'
seconds_check_out = 15.
care = 'no'
c_angle = 30

b_type = 'fgs'
e_folder = efs_folder

del_data, '*'

;;;; whether to translate to the propagate direction
v_suffix = '' ;;; X direction
;v_suffix = '_vperp'
;v_suffix = '_vperp2'
;v_suffix = '_vxy'
;v_suffix = '_vdht'
vc_use = 'vperp2'

;list_name = 'dfb_list_lead_tail_2013_tailward_df_binxbout.txt'
list_name = 'dfb_list_lead_tail.txt'

;;; load events
events = load_list(list_name, folder = listfolder)
;; diagnose
events = events(*, 100:110)

for i = 0, n_elements(events[0,*])-1 do begin
	time = time_double(events[0,i])
	sc = events[3,i]
	probe = sc
	trange = [time-180., time+180.]

	;;;;; test electric field data
	;thick = df_thick(events[*,i], normal_method = method, v_use = v_use, seconds_check_in = 15., seconds_check_out = seconds_check_out, care = care, b_type = b_type, fgs_folder = fgs_folder, fgl_folder = fgl_folder, e_folder = e_folder, /reduce, c_angle = c_angle, maxgap = 10)
	;copy_data, 'th'+sc+'_efs_dsl', 'th'+sc+'_efs_dsl_old'
	;del_data, 'th'+sc+'_efs_dsl'
;	;thm_load_fit, probe = sc, trange = [time-180., time+180.], coord='dsl', level = 2, data='efs'
	;thm_load_fit, probe = sc, trange = [time-180., time+180.], coord='dsl', level = 1, data='efs', suffix = '_dsl'
	;thm_load_fit, probe = sc, trange = [time-180., time+180.], coord='gsm', level = 1, data='fgs'
	;load_bin_data, probe = sc, trange = [time-180., time+180.], datatype = 'eff_dot0_dsl', datafolder = eff_folder
	;copy_data, 'th'+sc+'_eff_dot0_dsl', 'th'+sc+'_eff_dot0_dsl_old'
	;del_data, 'th'+sc+'_eff_dot0_dsl'
	;thm_load_efi, probe = sc, trange = [time-180., time+180.], level=2, coord='dsl', data='eff_dot0'
	;copy_data, 'th'+sc+'_eff_dot0_dsl', 'th'+sc+'_eff_dot0_dsl_l2'
	;thm_load_efi, probe = sc, trange = [time-180., time+180.], level=1, coord='dsl', data='eff_dot0'
	;copy_data, 'th'+sc+'_eff_dot0', 'th'+sc+'_eff_dot0_dsl_l1'
	;tplot, ['th'+sc+'_'+b_type+'_gsm', 'th'+sc+'_'+b_type, 'th'+sc+'_efs_dsl_old', 'th'+sc+'_efs_gsm', 'th'+sc+'_efs_dsl', 'th'+sc+'_eff_dot0_dsl_old', 'th'+sc+'_eff_dot0_dsl_l2', 'th'+sc+'_eff_dot0_dsl_l1'], trange = [time-180., time+180.]
	;timebar, time, line = 1

	;;;;; test Vi new vs old

	;;;; old
	load_bin_data, probe = probe, trange = trange, datatype = 'vi', datafolder = vi_folder
	options, 'th'+probe+'_ptix_velocity_gsm', colors = [2,4,6]

	;;;; new
	esa_datatype = 'peir'
	sst_datatype = 'psif'
	thm_load_state, probe = probe, trange = trange, /get_sup
	;;; if using automaticly determined bad bins, do not set sst_sun_bins
	combined = thm_part_combine(probe=probe, trange=trange, $
	                            esa_datatype=esa_datatype, sst_datatype=sst_datatype, $
	                            orig_esa=esa, orig_sst=sst)
	thm_part_products, dist_array=combined, outputs='moments'
	
	thm_cotrans, 'th'+probe+'_ptirf_velocity', 'th'+probe+'_ptirf_velocity_gsm', in_coord = 'dsl', out_coord = 'gsm'

	;ylim, 'th'+probe+'_ptirf_velocity_gsm', -300, 300
	;ylim, 'th'+probe+'_ptix_velocity_gsm', -300, 300
	
	time_stamp,/off
	
	tplot, ['th'+probe+'_ptix_velocity_gsm', 'th'+probe+'_ptirf_velocity_gsm'], trange = trange
	timebar_mass, [0, 100., 200., 300., 400., 600., 800., 1000], varname = ['th'+probe+'_ptix_velocity_gsm', 'th'+probe+'_ptirf_velocity_gsm'], /databar
	makepng, pic_folder+'/events_temp/event_'+strcompress(string(i), /remove)
endfor

stop
end
