pro main_multi
;; find out multi-point observations and do things
;; try to infer the expansion of the bundles with multi-point observations
thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

;mark = 'all'
mark = 'earth'

cc_req = 0.5 ;; for deHoffman teller method

secs_multi = 30. ;; for paper

year_suf = '' ;; 2007-2011 seasons
;year_suf = '_2013' ;; 2013 season

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

;;; time range to get t_out (minus this seconds)
t_out_suf = '' ;; default: 15s
;t_out_suf = '_0'
;t_out_suf = '_3'
;t_out_suf = '_5'

;;; set the range of events
xbin_bounds = [-30., -6.] ;; all
;xbin_bounds = [-15., -12., -11, -10, -9, -8, -6] ;; highest resolution, used for ratio_convex, 15s.
;xbin_bounds = [-15., -11., -8.5, -6] ;; for radius, 15s 
;xbin_bounds = [-15., -11., -10., -8.5, -6] ;; for dy, vperp2 as x seems to be better

;;; choose method for normal direction
;method_suffix = 'minvar'
method_suffix = 'binxbout'

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

;;;; which velocity to take when computing the expansion/contraction
;vc_use = 'vperp' ;; given by df_thick, experienced interpolation
;vc_use = 'vperp2' ;; exactly the vperp over the few points of the DF.
vc_use = 'vbh' ;; the vperp behind the DF.
;vc_use = 'vdht'

;;;; for vbh only: the seconds behind to count as behind the DF.
secs_bh_start = 0.
secs_bh_end = 7.

;;;; how to judge expand or contract
;judge_reshape_para = 'dvy'
;judge_reshape_para = 'angle'
judge_reshape_para = 'both'

;;;; whether to print details (dy, dn angles)
print_details = 'yes'
;print_details = 'no'

;;;; whether to translate to the propagate direction
v_suffix = '' ;;; X direction
;v_suffix = '_vperpX'
;v_suffix = '_vperp2X'
;v_suffix = '_vxyX'
;v_suffix = '_vdhtX'

;;;; max dx and dz in RE
;dx_max = 10. ;; max is 2.81 RE for 15s
;dz_max = 10. ;; max is 1.76RE for 15s
dx_max = 3.
dz_max = 2.

;;; revise seconds_check to test
seconds_check_in = 15.
case t_out_suf of
'': seconds_check_out = 15.
'_0': seconds_check_out = 0.
'_3': seconds_check_out = 3.
'_5': seconds_check_out = 5.
endcase

if strcmp(mark, 'all') then begin
	list_suffix = '_fgs'
endif else begin
	;;; select method to get the n directions
	case method_suffix of
		'minvar': method_num = 1
		'binxl': method_num = 2
		'binxbout': method_num = 3
	endcase
	;;; choose liberal dataset or not ;;;;;;;
	liberal_suf = '' ;; not liberal
	;liberal_suf = '_liberal' ;; liberal
	
	if strcmp(liberal_suf, '_liberal') then begin
		savefolder = savefolder+'/liberal'
		c_angle = 10.
	endif else begin
		c_angle = 30.
	endelse

	list_suffix = year_suf+'_'+mark+'ward_df_'+method_suffix+liberal_suf+t_out_suf

	;;;; for saving data
	if strcmp(two_side, 'yes') then suffix2 = '_two' else suffix2 = ''
	if strcmp(large_angle, 'yes') then suffix3 = '_largen' else suffix3 = ''
	if strcmp(large_dy, 'yes') then suffix4 = '_largey' else suffix4 = ''
endelse

;;; load events
events = load_list('dfb_list_lead_tail'+list_suffix+'.txt', folder = listfolder)
;; diagnose
;events = events(*, 0:3)

pos_list = event_status_instant(events, 'pos', time_length=3, datafolder=pos_folder)
pos = pos_list(2:4, *)
dataout_simple, save_folder+'/pos'+list_suffix, pos
;pos = datain_simple(save_folder+'/pos'+list_suffix+'.dat', dim=3, type='double')
x = pos(0,*)/RE
y = pos(1,*)/RE
z = pos(2,*)/RE

;; arrays for results
n_concave_arr = intarr(n_elements(xbin_bounds)-1)+!values.f_nan
n_convex_arr = n_concave_arr
n_expand_arr = n_convex_arr
n_contract_arr = n_convex_arr
r_mean_arr = fltarr(n_elements(xbin_bounds)-1)+!values.f_nan
r_std_arr = r_mean_arr
r_stat_arr = fltarr(n_elements(xbin_bounds)-1, 3)+!values.f_nan ;; lower quartile, median, higher quartile
dy_stat = r_stat_arr
dy_concave_stat = r_stat_arr
dy_convex_stat = r_stat_arr
dy_expand_stat = r_stat_arr
dy_contract_stat = r_stat_arr

t_out_cmpr_arr = [0d] ;; used to compare t_out from different methods
t_0_cmpr_arr = [0d] ;; to compare t_out with

;;; do a looop for x bins.
for i_dis = 0, n_elements(xbin_bounds)-2 do begin
	xrange = xbin_bounds(i_dis:i_dis+1)
	i_this = where((x gt xrange(0)) and (x lt xrange(1)), j_this)
	if j_this gt 1 then begin
		events_this = events(*, i_this)
		;;; get ny-y relation or probe distance distribution ;;;;;
		;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; get multi_point observations ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
		multi_list = multi_probe(events_this, c_seconds = secs_multi)
		
		;;;;; get the distances and Y positions ;;;;;
		multi_dxyz = [0., 0., 0.]
		multi_x = [0., 0.]
		multi_y = [0., 0.]
		multi_events = ['time', 'sc', 'time', 'sc']
		multi_Bx0 = [!values.f_nan, !values.f_nan]
		multi_By0 = [!values.f_nan, !values.f_nan]
		multi_vc_y = [!values.f_nan, !values.f_nan]
		multi_ny = [!values.f_nan, !values.f_nan]
		if ~strcmp(mark, 'all') then begin
			multi_nfront = [0., 0.]
			multi_np = [0., 0.]
			multi_dn = [0., 0., 0.]
			multi_vp = [0., 0., 0.]
			multi_vy = [0., 0.]
			multi_vx = [0., 0.]
			multi_dp = [0.] ;; location in the direction of average propagation
			multi_dfront = [0.] ;; location in the DF line
			multi_dnfront = [0.]
			multi_dvfront = [0.]
			multi_gamma = [0.] ;; the angle differences between the two azimuths (in XY)
			multi_r = [0.]
			multi_anglev = [0.]
		endif
		for i = 0, n_elements(multi_list(0,*))-1 do begin
			j_have = where(~strcmp(multi_list(*,i), 'no'), n_j_have)
			if n_j_have gt 0 then begin
				print, multi_list(j_have, i)
				times = time_double(multi_list(j_have, i))
				probes = probes_num_string(j_have)
				del_data, '*'
				pos = dblarr(3, n_elements(j_have))
				if ~strcmp(mark, 'all') then begin
					normal = dblarr(3, n_elements(j_have))
					bq = dblarr(3, n_elements(j_have))
					vp = dblarr(3, n_elements(j_have))
					v_convect = dblarr(3, n_elements(j_have))
				endif
				for j = 0, n_elements(j_have)-1 do begin
					;; get the position
					load_bin_data, probe = probes(j), trange = [times(j)-90., times(j)+90.], datatype='pos', /tclip, datafolder = pos_folder
					get_data, 'th'+probes(j)+'_state_pos_tclip', data = data
					pos(0,j) = mean(data.y(*,0), /nan)/RE
					pos(1,j) = mean(data.y(*,1), /nan)/RE
					pos(2,j) = mean(data.y(*,2), /nan)/RE
					if ~strcmp(mark, 'all') then begin
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
							if strcmp(vc_use, 'vbh') then begin
								v_no_meaning = vdht_calc([t_in+secs_bh_start, t_in+secs_bh_end], probe_arr = events_do(3,*), btype = btype, etype = etype, b_folder = b_folder, e_folder = e_folder, maxgap = 10, cc = cc, vperp = vbh)
								v_convect(*,j) = vbh
							endif else begin
								vdht = vdht_calc([t_out, t_in], probe_arr = events_do(3,*), btype = btype, etype = etype, b_folder = b_folder, e_folder = e_folder, maxgap = 10, cc = cc, vperp = vperp2)
								t_out_cmpr_arr = [t_out_cmpr_arr, t_out]
								t_0_cmpr_arr = [t_0_cmpr_arr, time_double(events_do(0))]
								if cc lt cc_req then vdht = !values.f_nan
								case vc_use of
								'vperp2': v_convect(*,j) = vperp2
								'vdht': v_convect(*,j) = vdht
								endcase
							endelse
						endelse
						case v_suffix of
						'': vp(*,j) = [100., 0., 0.]
						'_vperpX': vp(*,j) = vperp
						'_vperp2X': vp(*,j) = vperp2
						'_vxyX': vp(*,j) = vxy
						'_vdhtX': vp(*,j) = vdht
						endcase
						normal(*,j) = n_this
						bq(*,j) = bq_this
					endif ;; if for mark not "all"
				endfor ;; for of j, all dual-observations

				;;; calculate the distance
				for j = 0, n_elements(j_have)-2 do begin
					for k = j+1, n_elements(j_have)-1 do begin
						;; events
						multi_events = [[multi_events], [time_string(times(j)), probes(j), time_string(times(k)), probes(k)]]
						multi_x = [[multi_x], [pos(0,j), pos(0,k)]]
						multi_y = [[multi_y], [pos(1,j), pos(1,k)]]
						print, [probes(j), probes(k)]
						print, [pos(0,j), pos(0,k)]
		;				if strcmp(probes(j), 'b') and strcmp(probes(k), 'c') then stop
						;; dxyz
						dxyz_this = pos(*, j)-pos(*, k)
						if ~strcmp(mark, 'all') then begin
							;;; velocity
							vtube = [mean([vp(0,j), vp(0,k)]), mean([vp(1,j), vp(1,k)]), mean([vp(2,j), vp(2,k)])] ;; hear, vp means the vector
							multi_vp = [[multi_vp], [vtube]]
							vp_strgh_xy = sqrt(total(vtube(0:1)^2))
							vdir_xy = vtube(0:1)/vp_strgh_xy
							front_dir_xy = [-vdir_xy(1), vdir_xy(0)]
							;;; dxyz
							dp = dxyz_this(0)*vdir_xy(0)+dxyz_this(1)*vdir_xy(1)
							multi_dp = [multi_dp, dp]
							dfront = dxyz_this(0)*front_dir_xy(0)+dxyz_this(1)*front_dir_xy(1)
							multi_dfront = [multi_dfront, dfront]
							dxyz_this = [dp, dfront, dxyz_this(2)]
							;;; angle of v_convect from dusk event to dawn event
							vxy_j = [v_convect(0:1, j), 0.]
							vxy_k = [v_convect(0:1, k), 0.]
							anglev_j2k = angle_vectors(vxy_k, vxy_j, about = 'z') ;; angle from Vk to Vj, positive when counter-clockwise when viewed from Z=+inf
							dvfront = (vxy_j(0)-vxy_k(0))*front_dir_xy(0)+(vxy_j(1)-vxy_k(1))*front_dir_xy(1)
							multi_vx = [[multi_vx], [v_convect(0, j), v_convect(0, k)]]
							multi_vy = [[multi_vy], [v_convect(1, j), v_convect(1, k)]]
							;;; dn
							dn_this = normal(*,j)-normal(*,k)
							dnfront = dn_this(0)*front_dir_xy(0)+dn_this(1)*front_dir_xy(1)
							np_j = normal(0,j)*vdir_xy(0)+normal(1,j)*vdir_xy(1)
							np_k = normal(0,k)*vdir_xy(0)+normal(1,k)*vdir_xy(1)
							nfront_j = normal(0,j)*front_dir_xy(0)+normal(1,j)*front_dir_xy(1)
							nfront_k = normal(0,k)*front_dir_xy(0)+normal(1,k)*front_dir_xy(1)
							multi_nfront = [[multi_nfront], [nfront_j, nfront_k]]
							multi_np = [[multi_np], [np_j, np_k]]
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
									dvfront = !values.f_nan
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
							multi_dxyz = [[multi_dxyz], [dxyz_this]]
							multi_dn = [[multi_dn], [dn_this]]
							multi_dnfront = [multi_dnfront, dnfront]
							multi_gamma = [multi_gamma, gamma_jk]
							multi_r = [multi_r, r_this]
							multi_anglev = [multi_anglev, anglev_j2k]
							multi_dvfront = [multi_dvfront, dvfront]
							multi_ny = [[multi_ny], [nfront_j, nfront_k]]
							multi_vc_y = [[multi_vc_y], [vxy_j(1), vxy_k(1)]]
							multi_Bx0 = [[multi_Bx0], [bq(0,j), bq(0,k)]]
							multi_By0 = [[multi_By0], [bq(1,j), bq(1,k)]]
						endif ;; if for mark not "all"
					endfor ;; for of k, mutral
				endfor ;; for of j, all dual-observations
			endif ;; if for have dual-observations
		endfor ;; for of i, element on multi_list
		
		if n_elements(multi_dxyz(0,*)) gt 1 then begin
			multi_events = multi_events(*, 1:*)
			multi_dxyz = multi_dxyz(*, 1:*)
			multi_x = multi_x(*, 1:*)
			multi_y = multi_y(*, 1:*)
			if ~strcmp(mark, 'all') then begin
				multi_dn = multi_dn(*, 1:*)
				multi_vp = multi_vp(*, 1:*)
				multi_vx = multi_vx(*, 1:*)
				multi_vy = multi_vy(*, 1:*)
				multi_nfront = multi_nfront(*, 1:*)
				multi_np = multi_np(*, 1:*)
				multi_dp = multi_dp(1:*)
				multi_dfront = multi_dfront(1:*)
				multi_dnfront = multi_dnfront(1:*)
				multi_dvfront = multi_dvfront(1:*)
				multi_gamma = multi_gamma(1:*) ;; the difference between n directions
				multi_r = multi_r(1:*)
				multi_anglev = multi_anglev(1:*)
				multi_ny = multi_ny(*, 1:*)
				multi_vc_y = multi_vc_y(*, 1:*)
				multi_Bx0 = multi_Bx0(*, 1:*)
				multi_By0 = multi_By0(*, 1:*)
			endif
		endif

		;;; only when multi_dxyz and multi_dyodny can match each other
		within_dp = abs(multi_dp) le dx_max 
		within_dfront = (abs(multi_dfront) le 100.) and (abs(multi_dfront) ge -0.1)
		within_dz = abs(multi_dxyz(2,*)) le dz_max 
		i_within = where(within_dp and within_dfront and within_dz, j_within)
		
		if j_within gt 0 then begin
			multi_events = multi_events(*, i_within)
			multi_dxyz = multi_dxyz(*, i_within)
			multi_x = multi_x(*, i_within)
			multi_y = multi_y(*, i_within)
			if ~strcmp(mark, 'all') then begin
				multi_dn = multi_dn(*, i_within)
				multi_vp = multi_vp(*, i_within)
				multi_vx = multi_vx(*, i_within)
				multi_vy = multi_vy(*, i_within)
				multi_nfront = multi_nfront(*, i_within)
				multi_np = multi_np(*, i_within)
				multi_dfront = multi_dfront(i_within)
				multi_dnfront = multi_dnfront(i_within)
				multi_dvfront = multi_dvfront(i_within)
				multi_gamma = multi_gamma(i_within)
				multi_r = multi_r(i_within)
				multi_anglev = multi_anglev(i_within)
				multi_ny = multi_ny(*, i_within)
				multi_vc_y = multi_vc_y(*, i_within)
				multi_Bx0 = multi_Bx0(*, i_within)
				multi_By0 = multi_By0(*, i_within)
			endif
			
			;;; tell concave or convex: (y1-y2)/(ny1-ny2)
			if ~strcmp(mark, 'all') then begin
				multi_dyodny = multi_dfront/multi_dnfront
				multi_dnyody = multi_dnfront/multi_dfront
				i_convex = where(multi_dyodny gt 0, j_convex) ;; convex
				i_concave = where(multi_dyodny lt 0, j_concave) ;; concave

				;;;; concave events
				print, 'Number of concave dual-observations:'
				print, j_concave
				n_concave_arr(i_dis) = j_concave
				if j_concave gt 0 then begin
					print, '-- Median dny/dy for concave dual-observations:'
					print, median(multi_dnyody(i_concave), /even)
					dfront_concave = multi_dfront(i_concave)
					;;;;; dy and dn
					rstat, abs(dfront_concave), med_dfront_cave, hing1_dfront_cave, hing2_dfront_cave
					dy_concave_stat(i_dis, *) = [[hing1_dfront_cave], [med_dfront_cave], [hing2_dfront_cave]]
					gamma_concave = multi_gamma(i_concave)
					rstat, abs(gamma_concave), med_gamma_cave, hing1_gamma_cave, hing2_gamma_cave
					if strcmp(print_details, 'yes') then begin
						print, 'For all concave events:'
						print, '-- Meidan dn angle:'
						print, median(abs(gamma_concave), /even)*180./!pi
						print, '-- Upper quartile dn angle:'
						print, hing2_gamma_cave*180./!pi
						print, '-- Maximum dn angle:'
						print, max(gamma_concave, /nan)*180./!pi
						print, '-- Median dy (d front direction):'
						print, median(abs(dfront_concave), /even)
						print, '-- Maximum dy (d front direction):'
						print, max(abs(dfront_concave), /nan)
					endif

					;;;;;;;;; check other aspects
					ny_concave = multi_ny(*, i_concave)
					vc_y_concave  = multi_vc_y(*, i_concave)
					Bx0_concave  = multi_Bx0(*, i_concave)
					By0_concave  = multi_By0(*, i_concave)
					;;;;; check the Pitkanen [2013] model
					vybx_crit = Bx0_concave(0,*)*vc_y_concave(0,*)*Bx0_concave(1,*)*vc_y_concave(1,*) gt 0.
					vyinc_crit = (abs(Bx0_concave(0,*)/Bx0_concave(1,*))-1)*(abs(vc_y_concave(0,*)/vc_y_concave(1,*))-1) gt 0.
					i_pitkanen_noby = where(vybx_crit and vyinc_crit, j_pitkanen_noby)
					print, '-- '+strcompress(string(j_pitkanen_noby), /remove)+' are consistent with Pitkanen [2013] picture (not assuming By control)'
					;;; same by direction can be assumed for solar wind By direction
					sameby_crit = By0_concave(0,*)*By0_concave(1,*) gt 0
					i_sameby = where(sameby_crit, j_sameby)
					print, '-- of all concave, '+strcompress(string(j_sameby), /remove)+' have the same By sign and can be assumed solar wind By.'
					vcby_crit = (Bx0_concave(0,*)*vc_y_concave(0,*)*By0_concave(0,*) gt 0) and (Bx0_concave(1,*)*vc_y_concave(1,*)*By0_concave(1,*) gt 0)
					i_pitkanen_sameby = where(vybx_crit and vyinc_crit and sameby_crit and vcby_crit, j_pikanen_sameby)
					print, '------ of these, '+strcompress(string(j_pikanen_sameby), /remove)+' are consistent with Pitkanen [2013] picture (assuming By control).'
				endif ;; if of concave events

				;;;; convex events
				print, 'Number of convex dual-observations:'
				print, j_convex
				n_convex_arr(i_dis) = j_convex
				if j_convex gt 0 then begin
					print, '-- Median dny/dy for convex dual-observations:'
					print, median(multi_dnyody(i_convex), /even)
					r_convex = multi_r(i_convex)
					dfront_convex = multi_dfront(i_convex)
					dnfront_convex = multi_dnfront(i_convex)
					dvfront_convex = multi_dvfront(i_convex)
					gamma_convex = multi_gamma(i_convex)
					anglev_convex = multi_anglev(i_convex)

					;;;;; infer the dvy of expand/contract
					dvy_infer_convex = dvfront_convex/dnfront_convex*2

					;;;;; dy and dn
					rstat, abs(dfront_convex), med_dfront_vex, hing1_dfront_vex, hing2_dfront_vex
					dy_convex_stat(i_dis, *) = [[hing1_dfront_vex], [med_dfront_vex], [hing2_dfront_vex]]
					gamma_convex = multi_gamma(i_convex)
					if strcmp(print_details, 'yes') then begin
						print, 'For all convex events:'
						print, '-- Meidan dn angle:'
						print, median(abs(gamma_convex), /even)*180./!pi
						print, '-- Maximum dn angle:'
						print, max(gamma_convex, /nan)*180./!pi
						print, '-- Median dy (d front direction):'
						print, median(abs(dfront_convex), /even)
						print, '-- Maximum dy (d front direction):'
						print, max(abs(dfront_convex), /nan)
					endif

					;;;;; scale size
					print, 'The scale size of all convex DFBs:'
					print, '-- Mean radius:'
					print, mean(r_convex)
					r_mean_arr(i_dis) = mean(r_convex)
					print, '-- Standard deviation of the radius:'
					print, stddev(r_convex)
					r_std_arr(i_dis) = stddev(r_convex)
					rstat, r_convex, med_r, hing1, hing2
					print, '-- Median radius:'
					print, med_r
					r_stat_arr(i_dis, *) = [[hing1], [med_r], [hing2]]

					;;;;;;;;; check other aspects
					ny_convex = multi_ny(*, i_convex)
					vc_y_convex  = multi_vc_y(*, i_convex)
					Bx0_convex  = multi_Bx0(*, i_convex)
					By0_convex  = multi_By0(*, i_convex)
					;;;;; check the Pitkanen [2013] model
					vybx_crit = Bx0_convex(0,*)*vc_y_convex(0,*)*Bx0_convex(1,*)*vc_y_convex(1,*) gt 0.
					vyinc_crit = (abs(Bx0_convex(0,*)/Bx0_convex(1,*))-1)*(abs(vc_y_convex(0,*)/vc_y_convex(1,*))-1) gt 0.
					i_pitkanen_noby = where(vybx_crit and vyinc_crit, j_pitkanen_noby)
					print, '-- '+strcompress(string(j_pitkanen_noby), /remove)+' are consistent with Pitkanen [2013] picture (not assuming By control)'
					;;; same by direction can be assumed for solar wind By direction
					sameby_crit = By0_convex(0,*)*By0_convex(1,*) gt 0
					i_sameby = where(sameby_crit, j_sameby)
					print, '-- of all convex, '+strcompress(string(j_sameby), /remove)+' have the same By sign and can be assumed solar wind By.'
					vcby_crit = (Bx0_convex(0,*)*vc_y_convex(0,*)*By0_convex(0,*) gt 0) and (Bx0_convex(1,*)*vc_y_convex(1,*)*By0_convex(1,*) gt 0)
					i_pitkanen_sameby = where(vybx_crit and vyinc_crit and sameby_crit and vcby_crit, j_pikanen_sameby)
					print, '------ of these, '+strcompress(string(j_pikanen_sameby), /remove)+' are consistent with Pitkanen [2013] picture (assuming By control).'

					;;;;;;; contract or expand
					if strcmp(judge_reshape_para, 'angle') then begin
						judge_expand = anglev_convex/dfront_convex gt 0
						judge_contra = anglev_convex/dfront_convex lt 0
					endif
					if strcmp(judge_reshape_para, 'dvy') then begin
						judge_expand = dvfront_convex/dfront_convex gt 0
						judge_contra = dvfront_convex/dfront_convex lt 0
					endif
					if strcmp(judge_reshape_para, 'both') then begin
						judge_expand = (dvfront_convex/dfront_convex gt 0) and (anglev_convex/dfront_convex gt 0)
						judge_contra = (dvfront_convex/dfront_convex lt 0) and (anglev_convex/dfront_convex lt 0)
					endif

					;;; expand
					i_expand = where(judge_expand, n_expand)
					print, 'Number of expanding events:'
					print, n_expand
					n_expand_arr(i_dis) = n_expand
					if n_expand gt 0 then begin
						dfront_expand = dfront_convex(i_expand)
						rstat, abs(dfront_expand), med_dfront_exp, hing1_dfront_exp, hing2_dfront_exp
						dy_expand_stat(i_dis, *) = [[hing1_dfront_exp], [med_dfront_exp], [hing2_dfront_exp]]
						gamma_expand = gamma_convex(i_expand)
						rstat, abs(gamma_expand), med_gamma_exp, hing1_gamma_exp, hing2_gamma_exp
						if strcmp(print_details, 'yes') then begin
							print, 'For all expand events:'
							print, '-- Meidan dn angle:'
							print, median(abs(gamma_expand), /even)*180./!pi
							print, '-- Upper quartile dn angle:'
							print, hing2_gamma_exp*180./!pi
							print, '-- Maximum dn angle:'
							print, max(gamma_expand, /nan)*180./!pi
							print, '-- Median dy (d front direction):'
							print, median(abs(dfront_expand), /even)
							print, '-- Maximum dy (d front direction):'
							print, max(abs(dfront_expand), /nan)
							if strcmp(judge_reshape_para, 'both') then begin
								rstat, dvy_infer_convex(i_expand), med_dvyinf_exp, hing1_dvyinf_exp, hing2_dvyinf_exp
								print, '-- Inferred expansion speed dVy:'
								print, 'lq:'+strcompress(string(hing1_dvyinf_exp))+'; med:'+strcompress(string(med_dvyinf_exp))+'; hq:'+strcompress(string(hing2_dvyinf_exp))
							endif
						endif

						;;;;;;;;; check other aspects
						ny_expand = ny_convex(*, i_expand)
						vc_y_expand = vc_y_convex(*, i_expand)
						Bx0_expand = Bx0_convex(*, i_expand)
						By0_expand = By0_convex(*, i_expand)
						;;;;; check the Pitkanen [2013] model
						vybx_crit = Bx0_expand(0,*)*vc_y_expand(0,*)*Bx0_expand(1,*)*vc_y_expand(1,*) gt 0.
						vyinc_crit = (abs(Bx0_expand(0,*)/Bx0_expand(1,*))-1)*(abs(vc_y_expand(0,*)/vc_y_expand(1,*))-1) gt 0.
						i_pitkanen_noby = where(vybx_crit and vyinc_crit, j_pitkanen_noby)
						print, '-- '+strcompress(string(j_pitkanen_noby), /remove)+' are consistent with Pitkanen [2013] picture (not assuming By control)'
						;;; same by direction can be assumed for solar wind By direction
						sameby_crit = By0_expand(0,*)*By0_expand(1,*) gt 0
						i_sameby = where(sameby_crit, j_sameby)
						print, '-- of all expand, '+strcompress(string(j_sameby), /remove)+' have the same By sign and can be assumed solar wind By.'
						vcby_crit = (Bx0_expand(0,*)*vc_y_expand(0,*)*By0_expand(0,*) gt 0) and (Bx0_expand(1,*)*vc_y_expand(1,*)*By0_expand(1,*) gt 0)
						i_pitkanen_sameby = where(vybx_crit and vyinc_crit and sameby_crit and vcby_crit, j_pikanen_sameby)
						print, '------ of these, '+strcompress(string(j_pikanen_sameby), /remove)+' are consistent with Pitkanen [2013] picture (assuming By control).'
					endif

					;;; contract
					i_contra = where(judge_contra, n_contra)
					print, 'number of contracting events:'
					print, n_contra
					n_contract_arr(i_dis) = n_contra
					if n_contra gt 0 then begin
						dfront_contra = dfront_convex(i_contra)
						rstat, abs(dfront_contra), med_dfront_con, hing1_dfront_con, hing2_dfront_con
						dy_contract_stat(i_dis, *) = [[hing1_dfront_con], [med_dfront_con], [hing2_dfront_con]]
						gamma_contra = gamma_convex(i_contra)
						rstat, abs(gamma_contra), med_gamma_con, hing1_gamma_con, hing2_gamma_con
						if strcmp(print_details, 'yes') then begin
							print, 'For all contract events:'
							print, '-- Meidan dn angle:'
							print, median(abs(gamma_contra), /even)*180./!pi
							print, '-- Upper quartile dn angle:'
							print, hing2_gamma_con*180./!pi
							print, '-- Maximum dn angle:'
							print, max(gamma_contra, /nan)*180./!pi
							print, '-- Median dy (d front direction):'
							print, median(abs(dfront_contra), /even)
							print, '-- Maximum dy (d front direction):'
							print, max(abs(dfront_contra), /nan)
							if strcmp(judge_reshape_para, 'both') then begin
								rstat, dvy_infer_convex(i_contra), med_dvyinf_con, hing1_dvyinf_con, hing2_dvyinf_con
								print, '-- Inferred contraction speed dVy:'
								print, 'lq:'+strcompress(string(hing1_dvyinf_con))+'; med:'+strcompress(string(med_dvyinf_con))+'; hq:'+strcompress(string(hing2_dvyinf_con))
							endif
						endif

						;;;;;;;;; check other aspects
						ny_contra = ny_convex(*, i_contra)
						vc_y_contra = vc_y_convex(*, i_contra)
						Bx0_contra = Bx0_convex(*, i_contra)
						By0_contra = By0_convex(*, i_contra)
						;;;;; check the Pitkanen [2013] model
						vybx_crit = Bx0_contra(0,*)*vc_y_contra(0,*)*Bx0_contra(1,*)*vc_y_contra(1,*) gt 0.
						vyinc_crit = (abs(Bx0_contra(0,*)/Bx0_contra(1,*))-1)*(abs(vc_y_contra(0,*)/vc_y_contra(1,*))-1) gt 0.
						i_pitkanen_noby = where(vybx_crit and vyinc_crit, j_pitkanen_noby)
						print, '-- '+strcompress(string(j_pitkanen_noby), /remove)+' are consistent with Pitkanen [2013] picture (not assuming By control)'
						;;; same by direction can be assumed for solar wind By direction
						sameby_crit = By0_contra(0,*)*By0_contra(1,*) gt 0
						i_sameby = where(sameby_crit, j_sameby)
						print, '-- of all contract, '+strcompress(string(j_sameby), /remove)+' have the same By sign and can be assumed solar wind By.'
						vcby_crit = (Bx0_contra(0,*)*vc_y_contra(0,*)*By0_contra(0,*) gt 0) and (Bx0_contra(1,*)*vc_y_contra(1,*)*By0_contra(1,*) gt 0)
						i_pitkanen_sameby = where(vybx_crit and vyinc_crit and sameby_crit and vcby_crit, j_pikanen_sameby)
						print, '------ of these, '+strcompress(string(j_pikanen_sameby), /remove)+' are consistent with Pitkanen [2013] picture (assuming By control).'
					endif ;; if of existing contract events
				endif else begin
					n_expand_arr(i_dis) = 0
					n_contract_arr(i_dis) = 0
				endelse ;; if/else of convex events
			endif ;; if of whether mark is "all"
		endif else begin
			print, 'MAIN_MULTI: No dual observation within the dX and dZ range!'
		endelse ;; if/else of whether within specified dX and dZ range
	endif else begin
		print, 'MAIN_MULTI: No event in this range!'
	endelse ;; if/else of whether in the x_bin range
endfor


if n_elements(xbin_bounds) le 2 then begin
	;;;;;; output data
	if strcmp(two_side, 'no') and strcmp(large_angle, 'no') and strcmp(large_dy, 'no') then begin
;		dataout_simple, save_folder+'/dnangle'+list_suffix, transpose(multi_gamma)
;		dataout_simple, save_folder+'/dfront'+list_suffix+v_suffix, transpose(multi_dfront)
;		dataout_simple, save_folder+'/dnfront'+list_suffix+v_suffix, transpose(multi_dnfront)
;		dataout_simple, save_folder+'/dxyz'+list_suffix, multi_dxyz
;		dataout_simple, save_folder+'/x'+list_suffix, multi_x
;		dataout_simple, save_folder+'/y'+list_suffix, multi_y
;		dataout_simple, save_folder+'/nfront'+list_suffix, multi_nfront
;		dataout_simple, save_folder+'/np'+list_suffix, multi_np
;		dataout_simple, save_folder+'/vx'+list_suffix, multi_vx
;		dataout_simple, save_folder+'/vy'+list_suffix, multi_vy
;		dataout_simple, save_folder+'/vp'+list_suffix, multi_vp
;		dataout_simple, save_folder+'/i_convex'+list_suffix+v_suffix, transpose(i_convex)
;		dataout_simple, save_folder+'/i_concave'+list_suffix+v_suffix, transpose(i_concave)
;		dataout_simple, save_folder+'/r_convex'+list_suffix+v_suffix, transpose(r_convex)
;		dataout_simple, save_folder+'/dvfront_convex'+list_suffix+v_suffix, transpose(dvfront_convex)
;		dataout_simple, save_folder+'/i_expand'+list_suffix+v_suffix+'_'+vc_use+'_judge_'+judge_reshape_para, transpose(i_expand) ;; of convex
;		dataout_simple, save_folder+'/i_contra'+list_suffix+v_suffix+'_'+vc_use+'_judge_'+judge_reshape_para, transpose(i_contra) ;; of convex
		stop
	endif
	;if n_elements(t_out_cmpr_arr) gt 1 then begin
	;	t_out_cmpr_arr = t_out_cmpr_arr(1:*)
	;	t_0_cmpr_arr = t_0_cmpr_arr(1:*)
	;	dataout_simple, save_folder+'/tout_dual'+list_suffix, transpose(t_out_cmpr_arr)
	;	dataout_simple, save_folder+'/t0_dual'+list_suffix, transpose(t_0_cmpr_arr)
	;endif
endif else begin
	;;;;;; plot the variations
	xbincntrs = 0.5*(xbin_bounds(1:*)+xbin_bounds(0:-2))
	;; convex ratio
	plot, xbincntrs, double(n_convex_arr)/double(n_convex_arr+n_concave_arr), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'Ratio: Convex/All'
	oplot, xbincntrs, double(n_convex_arr)/double(n_convex_arr+n_concave_arr), psym = 4
	makepng, pic_folder+'/ratio_convex'
	;; convex ratio
	plot, xbincntrs, double(n_expand_arr)/double(n_expand_arr+n_contract_arr), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'Ratio: Expand/Convex'
	oplot, xbincntrs, double(n_expand_arr)/double(n_expand_arr+n_contract_arr), psym = 4
	makepng, pic_folder+'/ratio_expand'
	;; radius of DFB
	plot, xbincntrs, r_stat_arr(*,1), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'DFB radius [RE]', color = 6
	oplot, xbincntrs, r_stat_arr(*,1), psym = 4, color = 6
	oplot, xbincntrs, r_stat_arr(*,0), color = 2
	oplot, xbincntrs, r_stat_arr(*,0), psym = 4, color = 2
	oplot, xbincntrs, r_stat_arr(*,2), color = 2
	oplot, xbincntrs, r_stat_arr(*,2), psym = 4, color = 2
	makepng, pic_folder+'/r_dfb'
	;; dy of convex/concave
	plot, xbincntrs, dy_convex_stat(*,1), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'dY [RE]', color = 6, title = 'Convex vs concave'
	oplot, xbincntrs, dy_convex_stat(*,1), psym = 4, color = 6
	oplot, xbincntrs, dy_convex_stat(*,0), color = 2
	oplot, xbincntrs, dy_convex_stat(*,0), psym = 4, color = 2
	oplot, xbincntrs, dy_convex_stat(*,2), color = 2
	oplot, xbincntrs, dy_convex_stat(*,2), psym = 4, color = 2
	oplot, xbincntrs, dy_concave_stat(*,1), line = 1, color = 6
	oplot, xbincntrs, dy_concave_stat(*,1), line = 1, psym = 4, color = 6
	oplot, xbincntrs, dy_concave_stat(*,0), line = 1, color = 2
	oplot, xbincntrs, dy_concave_stat(*,0), line = 1, psym = 4, color = 2
	oplot, xbincntrs, dy_concave_stat(*,2), line = 1, color = 2
	oplot, xbincntrs, dy_concave_stat(*,2), line = 1, psym = 4, color = 2
	makepng, pic_folder+'/dy_convex_concave'
	;; dy of expand/contract
	plot, xbincntrs,  dy_expand_stat(*,1), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'dY [RE]', color = 6, title = 'Expand vs contract'
	oplot, xbincntrs, dy_expand_stat(*,1), psym = 4, color = 6
	oplot, xbincntrs, dy_expand_stat(*,0), color = 2
	oplot, xbincntrs, dy_expand_stat(*,0), psym = 4, color = 2
	oplot, xbincntrs, dy_expand_stat(*,2), color = 2
	oplot, xbincntrs, dy_expand_stat(*,2), psym = 4, color = 2
	oplot, xbincntrs, dy_contract_stat(*,1), line = 1, color = 6
	oplot, xbincntrs, dy_contract_stat(*,1), line = 1, psym = 4, color = 6
	oplot, xbincntrs, dy_contract_stat(*,0), line = 1, color = 2
	oplot, xbincntrs, dy_contract_stat(*,0), line = 1, psym = 4, color = 2
	oplot, xbincntrs, dy_contract_stat(*,2), line = 1, color = 2
	oplot, xbincntrs, dy_contract_stat(*,2), line = 1, psym = 4, color = 2
	makepng, pic_folder+'/dy_expand_contract'
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop

end
