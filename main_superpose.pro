pro main_superpose
thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

;; for dawn and dusk
ny_crit = 0.2
x_sep = -12.

;; superpose time range
pretime = 3.
afttime = 3.

;; the type of normal direction 
;method_suffix = 'minvar'
method_suffix = 'binxbout'

;;;; which velocity to take when computing the expansion/contraction
;vc_use = 'vperp' ;; given by df_thick, experienced interpolation
vc_use = 'vperp2' ;; exactly the vperp over the few points of the DF.
;vc_use = 'vdht'

;;; magnetic field & electric field used for vdHT and vperp
btype = 'fgl'
b_folder = fgl_folder
;;; high-res, slower
;etype = 'eff'
;e_folder = eff_folder
;;; low-res, faster
etype = 'efs'
e_folder = efs_folder

;;; time range to get t_out (minus this seconds)
t_out_suf = '' ;; default: 15s
;t_out_suf = '_0'
;t_out_suf = '_3'
;t_out_suf = '_5'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; end of settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
case etype of
'eff': etype_tvname = 'eff_dot0'
else: etype_tvname = etype
endcase

list_suffix = '_earthward_df_'+method_suffix+t_out_suf
listname = 'dfb_list_lead_tail'+list_suffix
pos = datain_simple(save_folder+'/pos'+list_suffix+'.dat', dim=3, type='double')
n = datain_simple(save_folder+'/n'+list_suffix+'.dat', dim=3, type='double')

;;; load events
events = load_list(listname+'.txt', folder = listfolder)
;;; diagnose
;events = events(*, 100:130)
;pos = pos(*, 100:130)
;n = n(*, 100:130)

x = pos(0,*)/RE
ny = n[1,*]

near_x = x gt x_sep
far_x = x lt x_sep
evening = ny gt ny_crit
morning = ny lt -ny_crit

store_data, 'bz_near', data = {type: 'bz', i_good: where(near_x)}
store_data, 'bz_far', data = {type: 'bz', i_good: where(far_x)}
store_data, 'vy_near_evening', data = {type: 'vby', i_good: where(near_x and evening)}
store_data, 'vy_near_morning', data = {type: 'vby', i_good: where(near_x and morning)}
store_data, 'vy_far_evening',  data = {type: 'vby', i_good: where(far_x and evening)}
store_data, 'vy_far_morning',  data = {type: 'vby', i_good: where(far_x and morning)}

;sup_names = ['bz_near', 'bz_far', 'vy_near_evening', 'vy_near_morning', 'vy_far_evening', 'vy_far_morning']
;sup_names = ['vy_far_morning', 'vy_far_evening', 'vy_near_evening', 'vy_near_morning']
sup_names = ['vy_near_morning']

for i = 0, n_elements(sup_names)-1 do begin
	get_data, sup_names(i), data = sup
	case sup.type of
	'bz': data_name = 'fgl_gsm_z_superpose_quartiles'
	'vby': data_name = etype_tvname+'_vperp_gsm_y_superpose_quartiles'
	endcase
	del_data, data_name
	if sup.i_good(0) ne -1 then begin
		if strcmp(sup.type, 'bz') then begin
			superpose_data_fast, 'fgl', components='z', t0p_list=[events(0, sup.i_good), events(3,sup.i_good)], pre=pretime, after=afttime, datafolder=fgl_folder;, detrend_v=vy_dusk
			type_suf = '_fgl'
		endif
		if strcmp(sup.type, 'vby') then begin
			superpose_data_fast_efi, 'vperp_'+etype, components='y', t0p_list=[events(0, sup.i_good), events(3,sup.i_good)], pre=pretime, after=afttime, e_folder = e_folder, btype = btype, b_folder = b_folder, rtrange = [time_double(events(0, sup.i_good))-180., time_double(events(0, sup.i_good))-120.], /clean, maxgap = 10
			type_suf = '_'+etype
		endif
		dataout, data_name, filename = save_folder+'/superpose_'+sup_names(i)+type_suf+list_suffix
	endif
endfor

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
stop
end
