pro main_single
;; try to infer the expansion of the bundles with single-point observations
thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

;; earthward or tailward
mark = 'earthward'
;mark = 'tailward'

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
etype = 'eff'
e_folder = eff_folder
;;; low-res, faster
;etype = 'efs'
;e_folder = efs_folder

;;; time range to get t_out (minus this seconds)
t_out_suf = '' ;; default: 15s
;t_out_suf = '_0'
;t_out_suf = '_3'
;t_out_suf = '_5'

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; end of settings ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; select method to get the n directions
case method_suffix of
	'minvar': method_num = 1
	'binxl': method_num = 2
	'binxbout': method_num = 3
endcase

list_suffix = '_'+mark+'_df_'+method_suffix+t_out_suf
listname = 'dfb_list_lead_tail'+list_suffix

;;; load events
events = load_list(listname+'.txt', folder = listfolder)
;;; diagnose
;events = events(*, 100:130)

;;; event positions
;pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
;pos = pos_list(2:4, *)
;x = pos(0,*)/RE
;y = pos(1,*)/RE
;z = pos(2,*)/RE

;;; set to use a particular set of events
;i_use = where(x lt 0.) ;; all
;;i_use = where((x gt -12) and (y lt -2)) ;; near, dawn
;;i_use = where((x gt -12) and (y gt 2)) ;; near, dusk
;;i_use = where(x lt -15) ;; far
;events = events(*, i_use)

if strcmp(vc_use, 'vperp') then begin
	no_use = df_thick(events, normal_method = method_suffix, v_use = 'efs', seconds_check_in = seconds_check_in, seconds_check_out = seconds_check_out, b_type = 'fgl', fgs_folder = fgs_folder, fgl_folder = fgl_folder, e_folder = efs_folder, normal_dir = n, maxgap = 10, vperp = vperp, vxy = vxy)
	v_convect = vperp
endif else begin
	;; get the normal
	n = df_normal(events, method = method_num, datatype=btype, bfolder=b_folder, care='yes', c_angle = 30., seconds_check_in = seconds_check_in, seconds_check_out = seconds_check_out, t_in = t_in, t_out = t_out)
	;;;; check t_out
	;i_close = where(t_out - (time_double(events(0,*))-15.) lt 3., n_close)
	;stop, n_close

	dataout_simple, save_folder+'/n'+list_suffix, n
	if n_elements(t_in) eq 0 then begin
		b_in_list = event_status_instant(events, btype, time_length = seconds_check_in/60., pre_time = -0.5*seconds_check_in/60., datafolder=b_folder, /max, comp_v = 2, time_interest = t_max_arr)
		t_in = t_max_arr(1,*)
	endif
	if n_elements(t_out) eq 0 then begin
		b_out_list = event_status_instant(events, btype, vtrange = [time_double(events(0,*))-seconds_check_out,t_max_arr(1,*)], datafolder=b_folder, /min, comp_v = 2, time_interest = t_min_arr)
  		t_out = t_min_arr(1,*)
	endif
	vdht = vdht_calc([t_out, t_in], probe_arr = events(3,*), btype = btype, etype = etype, b_folder = b_folder, e_folder = e_folder, maxgap = 10, cc = cc, vperp = vperp2)
	case vc_use of
	'vperp2': v_convect = vperp2
	'vdht': v_convect = vdht
	endcase
endelse

dataout_simple, save_folder+'/'+vc_use+'_'+etype+list_suffix, v_convect
;dataout_simple, save_folder+'/n'+list_suffix, n

;;;; components
;nx = n(0,*)
;ny = n(1,*)
;nz = n(2,*)
;vc_x = v_convect(0,*)
;vc_y = v_convect(1,*)
;vc_z = v_convect(2,*)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; modify for binned values ;;;;;;;;;;;;;;;;;;;
;;;;;; qtt 
;;qtt_bin = ny
;;qtt_title = 'ny'
;;qtt_str = 'ny'
;;vertical = 0
;;binsize = 0.25
;;qtt_range = [-1., 1.]
;;bin_boundaries = [-1., -0.5, 0., 0.5, 1.]
;
;qtt_bin = x
;qtt_title = 'X [RE]'
;qtt_str = 'x'
;vertical = 0
;binsize = 1.5
;qtt_range = [-6., -30.]
;;bin_boundaries = [-1., -0.5, 0., 0.5, 1.]
;
;;;;;; qtt 2
;;qtt_2_bin = vc_y
;;qtt_2_title = 'Vy [km/s]'
;;qtt_2_str = 'vy'
;;qtt_2_range = [-250., 250.]
;;k_c = 5
;
;qtt_2_bin = abs(vc_y)
;qtt_2_title = '|Vy| [km/s]'
;qtt_2_str = 'vyabs'
;qtt_2_range = [0., 300.]
;k_c = 5
;
;stat_plot, transpose(qtt_bin), transpose(qtt_2_bin), k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title, qtt_title = qtt_title, kinbin = kinbin, bincntrs_out = bincenters, pm = pm, ratio_pm = ratio_pm, vertical = vertical, bin_boundaries = bin_boundaries, title = title, avrg = avrg, std = std, med = med, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
;oplot, !x.crange, [0,0]
;oplot, [0,0], !y.crange
;makepng, pic_folder+'/'+qtt_2_str+'_'+qtt_str

stop
end
