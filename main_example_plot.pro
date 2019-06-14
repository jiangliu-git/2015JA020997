pro main_example_plot
;; plot an example of multi-point observation
thm_init
computer = 'I:'
;computer = '/home/jliu'
@folders

i_event_plot = 19 ;; for paper. this index is for large angle only, need to chagne large_angle to 'yes' and angle criterion only to get this.
;i_event_plot = 18 ;; for presentation

;plot_quantities = 'more'
plot_quantities = 'less'

;;; plot constants
if strcmp(plot_quantities, 'more') then begin
	xsize_tplot = 6
	ysize_tplot = 7.6
	abc = ['(a)','(b)','(c)','(d)','(e)','(f)']
	x_abc = 0.2
	y_abc = [0.9, 0.73, 0.58, 0.4, 0.29, 0.17]
	abc_pos = '(g)'
	;;;; right margin signs
	hpos = 0.9
	vert = 0.4
	hori = 0.03
	ypx_high = 0.75
	ypx_low = 0.3
endif
if strcmp(plot_quantities, 'less') then begin
	xsize_tplot = 6
	ysize_tplot = 6.6
	abc = ['(a)','(b)','(c)','(d)']
	x_abc = 0.2
	y_abc = [0.8, 0.67, 0.45, 0.23]
	abc_pos = '(e)'
	;;;; right margin signs
	hpos = 0.9
	vert = 0.4
	hori = 0.03
	ypx_high = 0.745
	ypx_low = 0.315
endif

xsize_pos = 4
ysize_pos = 3

;; some constants for loading data
minutes_load = 1.

;;;; whether to count only observations at two sides of the bubble
two_side = 'no'
;two_side = 'yes'

;;;; whether to count only observations having large theta difference
large_angle = 'no'
;large_angle = 'yes'
angle_req = 30.

;;;; whether to count only observations having large Y difference
large_dy = 'no'
;large_dy = 'yes'
dy_req = 0.5

;;; magnetic field & electric field used for vdHT and vperp
;;; high-res, slower
btype = 'fgl'
etype = 'eff'
b_folder = fgl_folder
e_folder = eff_folder
;;; low-res, faster
;btype = 'fgs'
;etype = 'efs'
;b_folder = fgs_folder
;e_folder = efs_folder

;;;; whether to translate to the propagate direction
v_suffix = '' ;;; X direction
;v_suffix = '_vperp'
;v_suffix = '_vperp2'
;v_suffix = '_vxy'
;v_suffix = '_vdht'
vc_use = 'vperp2'

;;; load events
events = load_list('dfb_list_lead_tail_earthward_df_binxbout.txt', folder = listfolder)
;; diagnose
;events = events(*, 0:3)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; get multi_point observations ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
multi_list = multi_probe(events, c_seconds = 30.)

expand_marker = 0

for i = 0, n_elements(multi_list(0,*))-1 do begin
	j_have = where(~strcmp(multi_list(*,i), 'no'), n_j_have)
	if n_j_have gt 0 then begin
		print, multi_list(j_have, i)
		times = time_double(multi_list(j_have, i))
		probes = probes_num_string(j_have)
		del_data, '*'
		pos = dblarr(3, n_elements(j_have))
		normal = dblarr(3, n_elements(j_have))
		t_in_arr = dblarr(n_elements(j_have))
		t_out_arr = dblarr(n_elements(j_have))
		vp = dblarr(3, n_elements(j_have))
		v_convect = dblarr(3, n_elements(j_have))
		for j = 0, n_elements(j_have)-1 do begin
			;; get the position
			load_bin_data, probe = probes(j), trange = [times(j)-90., times(j)+90.], datatype='pos', /tclip, datafolder = pos_folder
			get_data, 'th'+probes(j)+'_state_pos_tclip', data = data
			pos(0,j) = mean(data.y(*,0), /nan)/RE
			pos(1,j) = mean(data.y(*,1), /nan)/RE
			pos(2,j) = mean(data.y(*,2), /nan)/RE
			events_do = [time_string(times(j)), time_string(times(j)), 'm', probes(j)]
			if strcmp(vc_use, 'vperp') then begin
				no_use = df_thick(events_do, normal_method = method_suffix, v_use = etype, seconds_check = 15., b_type = btype, fgs_folder = fgs_folder, fgl_folder = fgl_folder, e_folder = e_folder, normal_dir = n_this, maxgap = 10, vperp = vperp, vxy = vxy)
				v_convect(*,j) = vperp
			endif else begin
				;; get the normal
				n_this = df_normal(events_do, method = method_num, datatype=btype, bfolder=b_folder, care='yes', c_angle = 30., seconds_check_in = seconds_check_in, seconds_check_out = seconds_check_out, t_in = t_in, t_out = t_out)
				bq_this_list = event_status_instant(events_do, 'fgs', time_length = 1., pre_time = 2.5, datafolder=fgs_folder)
				bq_this = bq_this_list(2:4,*)

				if n_elements(t_in) eq 0 then begin
					b_in_list = event_status_instant(events_do, btype, time_length = seconds_check_in/60., pre_time = -0.5*seconds_check_in/60., datafolder=b_folder, /max, comp_v = 2, time_interest = t_max_arr)
					t_in = t_max_arr(1,*)
				endif
				if n_elements(t_out) eq 0 then begin
					b_out_list = event_status_instant(events_do, btype, vtrange = [time_double(events_do(0,*))-seconds_check_out,t_max_arr(1,*)], datafolder=b_folder, /min, comp_v = 2, time_interest = t_min_arr)
					t_out = t_min_arr(1,*)
				endif
				vdht = vdht_calc([t_out, t_in], probe_arr = events_do(3,*), btype = btype, etype = etype, b_folder = b_folder, e_folder = e_folder, maxgap = 10, cc = cc, vperp = vperp2)
				case vc_use of
				'vperp2': v_convect(*,j) = vperp2
				'vdht': v_convect(*,j) = vdht
				endcase
			endelse
			case v_suffix of
			'': vp(*,j) = [100., 0., 0.]
			'_vperp': vp(*,j) = vperp
			'_vperp2': vp(*,j) = vperp2
			'_vxy': vp(*,j) = vxy
			'_vdht': vp(*,j) = vdht
			endcase
			normal(*,j) = n_this
			t_in_arr(j) = t_in
			t_out_arr(j) = t_out
		endfor ;; for of j, all dual-observations

		;;; calculate the distance
		for j = 0, n_elements(j_have)-2 do begin
			for k = j+1, n_elements(j_have)-1 do begin
				pos_j = pos(*, j)
				pos_k = pos(*, k)
				x_j = pos_j(0)
				y_j = pos_j(1)
				x_k = pos_k(0)
				y_k = pos_k(1)
				;; dxyz
				dxyz_this = pos_j-pos_k
				;;; velocity
				vp_ave = [mean([vp(0,j), vp(0,k)]), mean([vp(1,j), vp(1,k)]), mean([vp(2,j), vp(2,k)])]
				vp_strgh_xy = sqrt(total(vp_ave(0:1)^2))
				vdir_xy = vp_ave(0:1)/vp_strgh_xy
				front_dir_xy = [-vdir_xy(1), vdir_xy(0)]
				;;; dxyz
				dp = dxyz_this(0)*vdir_xy(0)+dxyz_this(1)*vdir_xy(1)
				dfront = dxyz_this(0)*front_dir_xy(0)+dxyz_this(1)*front_dir_xy(1)
				dxyz_this = [dp, dfront, dxyz_this(2)]
				;;; angle of v_convect from dusk event to dawn event
				vxy_j = [v_convect(0:1, j), 0.]
				vxy_k = [v_convect(0:1, k), 0.]
				anglev_j2k = angle_vectors(vxy_k, vxy_j, about = 'z') ;; angle from Vk to Vj, positive when counter-clockwise when viewed from Z=+inf
				dvfront = (vxy_j(0)-vxy_k(0))*front_dir_xy(0)+(vxy_j(1)-vxy_k(1))*front_dir_xy(1)
				;;; dn
				normal_j = normal(*,j)
				normal_k = normal(*,k)
				dn_this = normal(*,j)-normal(*,k)
				dnfront = dn_this(0)*front_dir_xy(0)+dn_this(1)*front_dir_xy(1)
				np_j = normal(0,j)*vdir_xy(0)+normal(1,j)*vdir_xy(1)
				np_k = normal(0,k)*vdir_xy(0)+normal(1,k)*vdir_xy(1)
				nfront_j = normal(0,j)*front_dir_xy(0)+normal(1,j)*front_dir_xy(1)
				nfront_k = normal(0,k)*front_dir_xy(0)+normal(1,k)*front_dir_xy(1)
				;;; the scale size
				alpha_j = atan2(nfront_j, np_j)
				alpha_k = atan2(nfront_k, np_k)
				gamma_jk = alpha_k-alpha_j
				;; get r
				;r_this = abs(dfront)/sqrt(2*(1-cos(gamma_jk))-(sin(alpha_j)-sin(alpha_k))^2)
				r_this = abs(dfront/(sin(alpha_j)-sin(alpha_k)))
				if strcmp(two_side, 'yes') then begin
					;;;; must be on the different sides to count ny
					if (sign(nfront_j)*sign(nfront_k) gt 0) or (abs(nfront_j-nfront_k) lt 0.0) then begin
						dn_this(1) = !values.f_nan 
						dnfront = !values.f_nan
						r_this = !values.f_nan
					endif
				endif
				if strcmp(large_angle, 'yes') then begin
					;;;; difference between two alphas should be less than 10 deg
					if abs(alpha_j-alpha_k) lt angle_req/180.*!pi then begin
						dn_this(1) = !values.f_nan 
						dnfront = !values.f_nan
						r_this = !values.f_nan
					endif
				endif

				;; judge to be convex
				convex_this = dfront/dnfront gt 0.
				concave_this = dfront/dnfront lt 0.
				;expand_this = anglev_j2k/dfront gt 0. ;; 'angle' criterion only
				expand_this = (anglev_j2k/dfront gt 0.) and (dvfront/dfront gt 0) ;; both criteria
				;;; draw figure if convex and expand
				if convex_this and expand_this then begin
					expand_marker = expand_marker+1
					;if expand_marker eq i_event_plot then begin ;; for paper
					;if abs(dnfront) gt 1.5 then begin ;; for check event
					if 1 gt 1.5 then begin ;; never do
						time_start = 0.5*(times(j)+times(k))
						trange_show = [time_start-minutes_load*60., time_start+minutes_load*60.] ;;; for DFB
						trange_load = [trange_show(0)-60., trange_show(1)+60.] ;;; for DFB

						;;;;; field data
						;; B fgl
						load_bin_data, trange = trange_load, probe = [probes(j), probes(k)], /tclip, datatype = 'fgl', datafolder = fgl_folder 
						;;; vperp_efs (2 probes)
						;load_efi_data, trange = trange_load, probe = probes(j), datatype = 'efs', rtrange = [times(j)-180., times(j)-120.], /tclip, e_folder = efs_folder, b_folder = fgl_folder
						;load_efi_data, trange = trange_load, probe = probes(k), datatype = 'efs', rtrange = [times(k)-180., times(k)-120.], /tclip, e_folder = efs_folder, b_folder = fgl_folder
						;; vperp_eff (2 probes)
						load_efi_data, trange = trange_load, probe = probes(j), datatype = 'eff', btype = 'fgl', rtrange = [times(j)-180., times(j)-120.], /tclip, e_folder = eff_folder, b_folder = fgl_folder
						load_efi_data, trange = trange_load, probe = probes(k), datatype = 'eff', btype = 'fgl', rtrange = [times(k)-180., times(k)-120.], /tclip, e_folder = eff_folder, b_folder = fgl_folder

						get_data, 'th'+probes[j]+'_eff_dot0_vperp_gsm_tclip', t, data
						store_data, 'th'+probes[j]+'_eff_dot0_vperp_gsm_tclip_xy', data={x:t, y:data[*, 0:1]}
						get_data, 'th'+probes[k]+'_eff_dot0_vperp_gsm_tclip', t, data
						store_data, 'th'+probes[k]+'_eff_dot0_vperp_gsm_tclip_xy', data={x:t, y:data[*, 0:1]}
						;; mark variables
						options, 'th?_fgl_gsm_tclip', colors = [2,4,6], ytitle = 'B!dGSM!n', ysubtitle = '!c[nT]', labels = ['B!dx!n', 'B!dy!n', 'B!dz!n'], labflag = 1
						options, 'th?_ef*_vperp_gsm_tclip', ytitle = 'V!dExB!n', ysubtitle = '!c[km/s]', colors = [2,4,6], labels = ['V!dx!n', 'V!dy!n', 'V!dz!n'], labflag = 1
						options, 'th?_ef*_vperp_gsm_tclip_xy', ytitle = 'V!dExB!n', ysubtitle = '!c[km/s]', colors = [2,4], labels = ['V!dx!n', 'V!dy!n'], labflag = 1
						if expand_marker eq i_event_plot then begin
							ylim, 'thd_eff_dot0_vperp_gsm_tclip_xy', -399., 499.
						endif else begin
							ylim, 'th?_ef*_vperp_gsm_tclip_xy', -999., 999.
						endelse

						;;;;;; plasma parameters
						;;; use saved data
						;;;; ni
						;;load_bin_data, trange = trange_load, probe = [probes(j), probes(k)], /tclip, datatype = 'ni', datafolder = ni_folder
						;;; Pall
						;load_bin_data, trange = trange_load, probe = [probes(j), probes(k)], /tclip, datatype = 'Pall', datafolder = Pall_folder
						;split_vec, 'th?_Pall_tclip'	
						;;options, 'th'+sc+'_ptix_density_tclip', ytitle = 'n!di!n', ysubtitle = '!c[cm!U-3!n]'
						;options, 'th?_Pall_tclip', ytitle = 'P', ysubtitle = '!c[nPa]', colors = [2,4,0], labels = ['P!db!n', 'P!dth!n', 'P!dttl!n'], labflag = 1
						;options, 'th?_Pall_tclip_y', ytitle = 'P!dth!n', ysubtitle = '!c[nPa]'

						;;; use newest data
						;;;;; determine the sst bins to remove manually. If use automatic bad bins, comment this part
						;thm_part_load,probe=probe,trange=trange,datatype=sst_datatype
						;thm_part_products,probe=probe,datatype=sst_datatype,trange=trange, sst_sun_bins = -1
						;tplot, 'th'+probe+'_psif_eflux_phi'
						;tm = gettime(/c)
						;edit3dbins, thm_part_dist(probe=probe, type='psif',/sst_cal), psif_badbins,/log
						;
						;;;;; new method, should be better
						;;;; if using automaticly determined bad bins, do not set sst_sun_bins
						;combined = thm_part_combine(probe=probe, trange=trange, $
						;                            esa_datatype=esa_datatype, sst_datatype=sst_datatype, $
						;                            orig_esa=esa, orig_sst=sst, sst_sun_bins=psif_badbins) 
						;thm_part_products, dist_array=combined, outputs='moments'
						;stop
	
						options, '*', thick = l_thick

						;title = strcompress(string(expand_marker), /remove)+' th'+probes(j)+' th'+probes(k)+' '+time_string(time_start)
						title = 'A Dual Observation Example'
						;;; order of plot names 
						if strcmp(plot_quantities, 'more') then begin
							qtts_j = ['th'+probes(j)+'_fgl_gsm_tclip', 'th'+probes(j)+'_Pall_tclip_y', 'th'+probes(j)+'_eff_dot0_vperp_gsm_tclip_xy']
							qtts_k = ['th'+probes(k)+'_fgl_gsm_tclip', 'th'+probes(k)+'_Pall_tclip_y', 'th'+probes(k)+'_eff_dot0_vperp_gsm_tclip_xy']
						endif
						if strcmp(plot_quantities, 'less') then begin
							qtts_j = ['th'+probes(j)+'_fgl_gsm_tclip', 'th'+probes(j)+'_eff_dot0_vperp_gsm_tclip_xy']
							qtts_k = ['th'+probes(k)+'_fgl_gsm_tclip', 'th'+probes(k)+'_eff_dot0_vperp_gsm_tclip_xy']
						endif
						if t_in_arr(j) le t_in_arr(k) then begin
							tp_names = [qtts_j, qtts_k] 
							ypx_j = ypx_high
							ypx_k = ypx_low
						endif else begin 
							tp_names = [qtts_k, qtts_j] 
							ypx_j = ypx_low
							ypx_k = ypx_high
						endelse
						;;; begin plotting
						popen, pic_folder+'/dual_example'+'_'+strcompress(string(expand_marker), /remove)
						print_options,xsize=xsize_tplot, ysize=ysize_tplot ;; use this for single plot
						tplot, tp_names, trange = trange_show, title = title
						timebar_mass, 0, varname=['th'+probes(j)+'_fgl_gsm_tclip', 'th'+probes(j)+'_eff_dot0_vperp_gsm_tclip_xy', 'th'+probes(k)+'_fgl_gsm_tclip', 'th'+probes(k)+'_eff_dot0_vperp_gsm_tclip_xy'], /databar, line = 3
						timebar_mass, [t_in_arr(j), t_out_arr(j)], varname=['th'+probes(j)+'_fgl_gsm_tclip', 'th'+probes(j)+'_eff_dot0_vperp_gsm_tclip_xy'], line = 1
						timebar_mass, [t_in_arr(k), t_out_arr(k)], varname=['th'+probes(k)+'_fgl_gsm_tclip', 'th'+probes(k)+'_eff_dot0_vperp_gsm_tclip_xy'], line = 1
						if strcmp(plot_quantities, 'more') then begin
							timebar_mass, [t_in_arr(j), t_out_arr(j)], varname='th'+probes(j)+'_Pall_tclip_y', line = 1
							timebar_mass, [t_in_arr(k), t_out_arr(k)], varname='th'+probes(k)+'_Pall_tclip_y', line = 1
						endif
	
						;;;; write abcs
						xyouts, replicate(x_abc, n_elements(tp_names)), y_abc, abc, /normal, charsize = 1.2
						;xyouts, 0.2, 0.7, strcompress(string(r_this), /remove), /normal
						make_px, [hpos, ypx_j], ver=vert, hor=hori, charname='P'+thm_probe_color(probes(j), /number), lchar=0.05, orientation = -90, alignment = 0.5, color = thm_probe_color(probes(j))
						make_px, [hpos, ypx_k], ver=vert, hor=hori, charname='P'+thm_probe_color(probes(k), /number), lchar=0.05, orientation = -90, alignment = 0.5, color = thm_probe_color(probes(k))
						pclose
						;makepng, pic_folder+'/events/dual_example_'+strcompress(string(expand_marker), /remove)

						;;;; plot the position/n/v
						;; infer the origin of the two DFB semicircles
						Oj = [x_j-r_this*cos(alpha_j), y_j-r_this*sin(alpha_j)]
						Ok = [x_k-r_this*cos(alpha_k), y_k-r_this*sin(alpha_k)]
						theta = findgen(30)/29.*!pi-0.5*!pi
						x_Oj = Oj(0)+r_this*cos(theta)
						y_Oj = Oj(1)+r_this*sin(theta)
						x_Ok = Ok(0)+r_this*cos(theta)
						y_Ok = Ok(1)+r_this*sin(theta)
						;; the tangential directions
						half_l_tang = 0.25
						tang_j = [-normal_j[1], normal_j[0]]
						tang_j = tang_j/sqrt(total(tang_j^2))
						tang_k = [-normal_k[1], normal_k[0]]
						tang_k = tang_k/sqrt(total(tang_k^2))
						x_tang_j = [x_j-tang_j[0]*half_l_tang*r_this, x_j+tang_j[0]*half_l_tang*r_this]
						y_tang_j = [y_j-tang_j[1]*half_l_tang*r_this, y_j+tang_j[1]*half_l_tang*r_this]
						x_tang_k = [x_k-tang_k[0]*half_l_tang*r_this, x_k+tang_k[0]*half_l_tang*r_this]
						y_tang_k = [y_k-tang_k[1]*half_l_tang*r_this, y_k+tang_k[1]*half_l_tang*r_this]
						;; arrowheads
						v_max = max([sqrt(total(vxy_j^2)), sqrt(total(vxy_j^2))], /nan)
						length_n = r_this/3.
						length_v = r_this/(1.5*v_max)
						x_n_heads = [x_j+length_n*normal_j(0), x_k+length_n*normal_k(0)]
						y_n_heads = [y_j+length_n*normal_j(1), y_k+length_n*normal_k(1)]
						x_v_heads = [x_j+length_v*vxy_j(0), x_k+length_v*vxy_k(0)]
						y_v_heads = [y_j+length_v*vxy_j(1), y_k+length_v*vxy_k(1)]
						;; decide the range of the plot
						x_min = min([x_Oj, x_Ok, x_tang_j, x_tang_k, x_n_heads, x_v_heads], /nan)
						x_max = max([x_Oj, x_Ok, x_tang_j, x_tang_k, x_n_heads, x_v_heads], /nan)
						x_span = x_max-x_min
						y_min = min([y_Oj, y_Ok, y_tang_j, y_tang_k, y_n_heads, y_v_heads], /nan)
						y_max = max([y_Oj, y_Ok, y_tang_j, y_tang_k, y_n_heads, y_v_heads], /nan)
						y_span = y_max-y_min
						rate_range = 0.15
						;;;; start plotting
						popen, pic_folder+'/dual_example_pos'+'_'+strcompress(string(expand_marker), /remove)
						print_options,xsize=xsize_pos, ysize=ysize_pos;; use this for single plot
						;title = strcompress(string(expand_marker), /remove)+' th'+probes(j)+' th'+probes(k)+' '+time_string(time_start)
						title = 'Satellite Positions'
						;; plot the circles
						plot, x_Oj, y_Oj, line = 2, xtitle = 'X [R!dE!n]', ytitle = 'Y [R!dE!n]', /isotropic, title = title, xrange = [x_max+rate_range*x_span, x_min], yrange = [y_max+rate_range*y_span, y_min-rate_range*y_span], xstyle = 1, ystyle = 1
						oplot, x_Ok, y_Ok, line = 2
						;; plot probes
						draw_circle, x_j, y_j, 0.05*r_this, /stroke, /fill, color_fill = thm_probe_color(probes(j))
						draw_circle, x_k, y_k, 0.05*r_this, /stroke, /fill, color_fill = thm_probe_color(probes(k))
						xyouts, x_j+0.02*(!x.crange(1)-!x.crange(0)), y_j, 'P'+thm_probe_color(probes(j), /number), /data, charsize = 0.7, color = thm_probe_color(probes(j))
						xyouts, x_k+0.02*(!x.crange(1)-!x.crange(0)), y_k, 'P'+thm_probe_color(probes(k), /number), /data, charsize = 0.7, color = thm_probe_color(probes(k))
						;; the tangential directions
						oplot, x_tang_j, y_tang_j
						oplot, x_tang_k, y_tang_k
						;; plot the arrows
						arrow, [x_j, x_k], [y_j, y_k], x_v_heads, y_v_heads, hsize = -0.3, /data, color = 85, thick = 2
						arrow, [x_j, x_k], [y_j, y_k], x_n_heads, y_n_heads, hsize = -0.3, /data
						;; make the unit vector
						arrow_legend, 0.1, 0.15, length = length_v, value = -100., v_str = '100', unit_str = 'km/s', qtt_str = 'V!dDF,XY!n', color = 85, y_up = 0.02, y_down = 0.06
						xyouts, 0.07*(!x.crange[1]-!x.crange[0])+!x.crange[0], 0.85*(!y.crange[1]-!y.crange[0])+!y.crange[0], abc_pos
						pclose
						;makepng, pic_folder+'/events/dual_example_'+strcompress(string(expand_marker), /remove)+'pos'
						
						;;;;;;;;; print velocity values ;;;
						;; infer the expansion speed
						;nxy_j = normal_j[0:1]
						;nxy_j = nxy_j/sqrt(total(nxy_j^2))
						;nxy_k = normal_k[0:1]
						;nxy_k = nxy_k/sqrt(total(nxy_k^2))
						;dvy_infer = (vxy_j[1]-vxy_k[1])/(nxy_j[1]-nxy_k[1])*2
						dvy_infer = (vxy_j[1]-vxy_k[1])/(normal_j[1]-normal_k[1])*2
						print, 'Vx ave: '+string(mean([vxy_j[0], vxy_k[0]]))
						print, 'dVy_infer: '+string(dvy_infer)
						print, 'Marker:'+string(expand_marker)
						print, 'ny 0:'+string(normal_j[1])
						print, 'nx 0:'+string(normal_j[0])
						print, 'ny 1:'+string(normal_k[1])
						print, 'nx 1:'+string(normal_k[0])
						;stop
					endif
				endif
			endfor ;; for of k, mutral
		endfor ;; for of j, all dual-observations
	endif ;; if for have dual-observations
endfor ;; for of i, element on multi_list

print, 'Marker final:'+string(expand_marker)

stop
end
