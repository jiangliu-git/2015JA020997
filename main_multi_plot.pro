pro main_multi_plot
;; plot things from main_multi
thm_init

computer = 'D:'
;computer = '/home/jliu'

@folders

;;; choose method for normal direction
;method_suffix = 'minvar'
method_suffix = 'binxbout'

;;;; whether to translate to the propagate direction
v_suffix = '' ;;; X direction
;v_suffix = '_vperpX'
;v_suffix = '_vperp2X'
;v_suffix = '_vxyX'
;v_suffix = '_vdhtX'

;;;; which velocity to take when computing the expansion/contraction
;vc_use = 'vperp' ;; given by df_thick, experienced interpolation
vc_use = 'vperp2' ;; exactly the vperp over the few points of the DF.
;vc_use = 'vdht'

;;;; decide whether use events with large dy only
;large_dy = 'yes'
large_dy = 'no'

;;;; how to judge expand or contract
;judge_reshape_para = 'dvy'
;judge_reshape_para = 'angle'
judge_reshape_para = 'both'

;;; time range to get t_out (minus this seconds)
t_out_suf = '' ;; default: 15s
;t_out_suf = '_0'

list_suffix = '_earthward_df_'+method_suffix+t_out_suf

x = datain_simple(save_folder+'/x'+list_suffix+'.dat', dim = 2, type = 'double')
y = datain_simple(save_folder+'/y'+list_suffix+'.dat', dim = 2, type = 'double')
nfront = datain_simple(save_folder+'/nfront'+list_suffix+'.dat', dim = 2, type = 'double')
np = datain_simple(save_folder+'/np'+list_suffix+'.dat', dim = 2, type = 'double')
vx = datain_simple(save_folder+'/vx'+list_suffix+'.dat', dim = 2, type = 'double')
vy = datain_simple(save_folder+'/vy'+list_suffix+'.dat', dim = 2, type = 'double')
dxyz = datain_simple(save_folder+'/dxyz'+list_suffix+'.dat',  dim = 3, type = 'double')
dnangle = transpose(datain_simple(save_folder+'/dnangle'+list_suffix+'.dat', dim = 1, type = 'double'))*180./!pi
dfront = transpose(datain_simple(save_folder+'/dfront'+list_suffix+v_suffix+'.dat',  dim = 1, type = 'double'))
dnfront = transpose(datain_simple(save_folder+'/dnfront'+list_suffix+v_suffix+'.dat',  dim = 1, type = 'double'))
r_convex = transpose(datain_simple(save_folder+'/r_convex'+list_suffix+v_suffix+'.dat', dim = 1, type = 'double'))
dvfront_convex = transpose(datain_simple(save_folder+'/dvfront_convex'+list_suffix+v_suffix+'.dat', dim = 1, type = 'double'))
i_convex = transpose(datain_simple(save_folder+'/i_convex'+list_suffix+v_suffix+'.dat', dim = 1, type = 'long'))
i_concave = transpose(datain_simple(save_folder+'/i_concave'+list_suffix+v_suffix+'.dat',dim = 1, type = 'long'))
i_expand = transpose(datain_simple(save_folder+'/i_expand'+list_suffix+v_suffix+'_'+vc_use+'_judge_'+judge_reshape_para+'.dat', dim = 1, type = 'long')) ;; of convex
i_contra = transpose(datain_simple(save_folder+'/i_contra'+list_suffix+v_suffix+'_'+vc_use+'_judge_'+judge_reshape_para+'.dat', dim = 1, type = 'long')) ;; of convex

nfront_xy = nfront/sqrt(nfront^2+np^2)
dnfront_xy = nfront_xy(1,*)-nfront_xy(0,*)
x_convex	= x(*, i_convex)
x_concave	= x(*, i_concave)
y_convex	= y(*, i_convex)
nfront_convex	= nfront(*, i_convex)
np_convex	= np(*, i_convex)
vy_convex	= vy(*, i_convex)
vx_convex	= vx(*, i_convex)
dfront_convex	= dfront(i_convex)
dfront_concave	= dfront(i_concave)
dnangle_convex	= dnangle(i_convex)
dnangle_concave	= dnangle(i_concave)
dnfront_convex	= dnfront(i_convex)
dnfront_xy_convex	= dnfront_xy(i_convex)
x_expand	= x_convex(*, i_expand)
x_contra	= x_convex(*, i_contra)
nfront_expand = nfront_convex(*, i_expand)
np_expand	= np_convex(*, i_expand)
vy_expand = vy_convex(*, i_expand)
vx_expand	= vx_convex(*, i_expand)
dfront_expand	= dfront_convex(i_expand)
dfront_contra	= dfront_convex(i_contra)
dnangle_expand	= dnangle_convex(i_expand)
dnangle_contra	= dnangle_convex(i_contra)
dvfront_expand	= dvfront_convex(i_expand)
dnfront_expand	= dnfront_convex(i_expand)
dnfront_xy_expand	= dnfront_xy_convex(i_expand)

help, where(abs(dnangle_expand) gt 26)
help, where(abs(dnfront_expand) gt 1.5)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; general quantities ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dxyz_abs = abs(dxyz)
;rstat, dxyz_abs(0,*), x_med, x_low, x_high
;rstat, dxyz_abs(1,*), y_med, y_low, y_high
;rstat, dxyz_abs(2,*), z_med, z_low, z_high
;print, 'All dual observations:'
;print, 'dX: '+'lowq='+strcompress(string(x_low),/remove)+' med='+strcompress(string(x_med),/remove)+' highq='+strcompress(string(x_high),/remove)
;print, 'dY: '+'lowq='+strcompress(string(y_low),/remove)+' med='+strcompress(string(y_med),/remove)+' highq='+strcompress(string(y_high),/remove)
;print, 'dZ: '+'lowq='+strcompress(string(z_low),/remove)+' med='+strcompress(string(z_med),/remove)+' highq='+strcompress(string(z_high),/remove)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; results with requirements ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dy_req = 0.5
;i_smalldy = where(abs(dfront) gt dy_req, n_smalldy)
;i_smalldy_convex = where(abs(dfront_convex) gt dy_req, n_smalldy_convex)
;print, n_smalldy
;print, n_smalldy_convex
;print, double(n_smalldy_convex)/double(n_smalldy)
;
;i_smalldy_expand = where(abs(dfront_expand) gt dy_req, n_smalldy_expand)
;i_smalldy_contra = where(abs(dfront_contra) gt dy_req, n_smalldy_contra)
;print, n_smalldy_expand
;print, n_smalldy_contra
;print, double(n_smalldy_expand)/double(n_smalldy_expand+n_smalldy_contra)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; dVy vs dny ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;dny_use = dnfront_expand
;ny_str = ''
;;dny_use = dnfront_xy_expand
;;ny_str = '!u*!n'
;
;if strcmp(t_out_suf, '') then begin
;	if strcmp(ny_str, '') then binbounds = [0., 0.1, 0.22, 0.37, 0.55, 1., 1.8] else binbounds = [0., 0.1, 0.3, 0.5, 1., 1.5, 2.]
;	if strcmp(ny_str, '') then ny_range = [0., 1.8] else ny_range = [0., 2.]
;	vy_range = [0., 299.]
;endif
;if strcmp(t_out_suf, '_0') then begin
;	binbounds = [0., 0.1, 0.32, 0.55,0.8, 1.7]
;	;binbounds = [0., 0.1, 0.36, 0.55, 1., 1.8]
;	;binbounds = [0., 0.25, 0.5, 1., 1.8]
;	ny_range = [0., 1.7]
;	vy_range = [0., 249.]
;endif
;
;title = 'Dual Observations of Expanding Events'
;
;popen, pic_folder+'/dual_dvy_dny'+t_out_suf
;print_options,xsize=3.8,ysize=3.
;stat_plot, abs(dny_use), abs(dvfront_expand), k_c = 3, bin_range = [0.,2.], binsize = 0.1, qtt_2_range = vy_range, qtt_range = ny_range, qtt_2_title = '|'+cap_delta_letter+'V!dDF,y!n| [km/s]', qtt_title = '|'+cap_delta_letter+'n!dy!n'+ny_str+'|', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick
;xyouts, !x.crange(0)+0.5*(!x.crange(1)-!x.crange(0)), !y.crange(1)+0.03*(!y.crange(1)-!y.crange(0)), title, align = 0.5, charsize = 0.9
;pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; Considering Slide, tube motion, and expansion ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;; infer the Vx of the flux tube assuming vy is completely from sliding motion: stat plot against dny
;npts_try = 1000
;max_try = 1500.
;;max_try = 200.
;vtubey_try = dindgen(npts_try)/(npts_try-1)*2*max_try-max_try
;i_check = where(abs(dnfront_expand) gt 1.)
;;i_check = where(abs(dnfront_expand) ge 0.) ;; all
;vtube_x_infer = dblarr(n_elements(i_check))
;vtube_y_infer = dblarr(n_elements(i_check))
;for i = 0, n_elements(i_check)-1 do begin
;	i_this = i_check[i]
;	;; probe 0
;	vys0 = vy_expand(0,i_this)-vtubey_try
;	vxs0 = -vys0*nfront_expand(0,i_this)/np_expand(0,i_this)
;	vtube_x0 = vx_expand(0,i_this)-vxs0
;	;; probe 1
;	vys1 = vy_expand(1,i_this)-vtubey_try
;	vxs1 = -vys1*nfront_expand(1,i_this)/np_expand(1,i_this)
;	vtube_x1 = vx_expand(1,i_this)-vxs1
;	;; get the max of the two
;	max_v = max([[vtube_x0], [vtube_x1]], dimension = 2)
;	;plot, vtubey_try, max_v, xtitle = 'Assumed V_tube_y [km/s]', ytitle = 'Required (higher) V_tube_x [km/s]', title = 'event'+strcompress(string(i_this))
;	vtube_x = min(max_v, i_min, /abs)
;	vtube_y = vtubey_try(i_min)
;	vtube_x_infer[i] = vtube_x
;	vtube_y_infer[i] = vtube_y
;	if (strcmp(t_out_suf, '') and (i eq 9)) or (strcmp(t_out_suf, '_0') and (abs(dnfront_expand[i_this]) gt 1.5)) then begin
;		print, 'Sc 0:'
;		print, 'Vobs,x:'+string(vx_expand(0,i))
;		print, 'Vobs,y:'+string(vy_expand(0,i))
;		print, 'Vs,x:'+string(vxs0[i_min])
;		print, 'Vs,y:'+string(vys0[i_min])
;		print, 'Vtube,x:'+string(vtube_x0[i_min])
;		print, 'Vtube,y:'+string(vtube_y)
;		print, 'Sc 1:'
;		print, 'Vobs,x:'+string(vx_expand(1,i))
;		print, 'Vobs,y:'+string(vy_expand(1,i))
;		print, 'Vs,x:'+string(vxs1[i_min])
;		print, 'Vs,y:'+string(vys1[i_min])
;		print, 'Vtube,x:'+string(vtube_x1[i_min])
;		print, 'Vtube,y:'+string(vtube_y)
;		stop
;	endif
;end
;
;;; plot values against dnfront
;ny_str = ''
;binbounds = [0., 0.1, 0.22, 0.37, 0.55, 1., 1.5, 2.]
;dnf_check = dnfront_expand(i_check)
;print, 'dny:'
;print, dnf_check
;print, 'Vtubex:'
;print, vtube_x_infer
;print, 'Vtubey:'
;print, vtube_y_infer
;stat_plot, abs(dnf_check), vtube_x_infer, k_c = 2, bin_range = [0.,2.], binsize = 0.1, qtt_2_range = vy_range, qtt_range = ny_range, qtt_2_title = 'V!dtube,x!n inferred [km/s]', qtt_title = '|'+cap_delta_letter+'n!dy!n'+ny_str+'|', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick
;;makepng, pic_folder+'/vtubex'+t_out_suf
;;;stat_plot, abs(dnf_check), vtube_y_infer, k_c = 3, bin_range = [0.,2.], binsize = 0.1, qtt_2_range = vy_range, qtt_range = ny_range, qtt_2_title = 'V!dtube,y!n inferred [km/s]', qtt_title = '|'+cap_delta_letter+'n!dy!n'+ny_str+'|', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick
;;;makepng, pic_folder+'/vtubey'
;;stat_plot, abs(dnf_check), abs(vtube_y_infer), k_c = 2, bin_range = [0.,2.], binsize = 0.1, qtt_2_range = vy_range, qtt_range = ny_range, qtt_2_title = '|V!dtube,y!n| inferred [km/s]', qtt_title = '|'+cap_delta_letter+'n!dy!n'+ny_str+'|', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick
;;makepng, pic_folder+'/vtubey_abs'+t_out_suf
;;;;;;;;;;;;;;;;
;
;;;;;;;;;;;;;;; what Vytube should be to have this speed for the tube Vx
;vtubex_should = 195. ; in km/s, for inferring expansion speed. This is from the small ny events
;;vtubex_should = 334. ; in km/s, for calculating vtubey_should. This is from the large dNy events
;;vtubex_should = 364. ; in km/s, for calculating vtubey_should. This is from the large dNy events
;
;vtubey_should_all = dblarr(n_elements(dnfront_expand))
;vexpand_should_all= dblarr(n_elements(dnfront_expand))
;for i = 0, n_elements(dnfront_expand)-1 do begin
;	;; probe 0
;	vxs0 = vx_expand(0,i)-vtubex_should
;	vys0 = -vxs0*np_expand(0,i)/nfront_expand(0,i)
;	vtubey_should0 = vy_expand(0,i)-vys0
;	;; probe 1
;	vxs1 = vx_expand(1,i)-vtubex_should
;	vys1 = -vxs1*np_expand(1,i)/nfront_expand(1,i)
;	vtubey_should1 = vy_expand(1,i)-vys1
;	;; the minimum one is for this
;	if abs(vtubey_should0) gt abs(vtubey_should1) then vtubey_should_this = vtubey_should1 else vtubey_should_this = vtubey_should0
;	vtubey_should_all(i) = vtubey_should_this
;	;; expansion speed
;	vexpand_should_all(i) = (vtubey_should0-vtubey_should1)/sign(dnfront_expand(i))
;
;	;;; examine
;	;print, 'sc0:'+string(vtubey_should0)
;	;print, 'sc1:'+string(vtubey_should1)
;
;	;; take
;	;if strcmp(vcomp, 'vx') then begin
;	;	oplot, nfront_expand[*, i], vx_expand[*, i]
;	;	oplot, nfront_expand[*, i], vx_expand[*, i], psym = 4
;	;endif
;	;if strcmp(vcomp, 'vy') then begin
;	;	oplot, nfront_expand[*, i], vy_expand[*, i]
;	;	oplot, nfront_expand[*, i], vx_expand[*, i], psym = 4
;	;endif
;	;if strcmp(vcomp, 'vtubex') and (n_elements(vtubey_try) lt 2) then begin
;	;	oplot, nfront_expand[*, i], [vtube_x0, vtube_x1]
;	;	oplot, nfront_expand[*, i], [vtube_x0, vtube_x1], psym = 4
;	;endif
;
;	;if abs(dnfront_expand(i)) gt 1.3 then stop
;	;if i eq 9 then begin
;	;	print, 'Sc 0:'
;	;	print, 'Vobs,x:'+string(vx_expand(0,i))
;	;	print, 'Vobs,y:'+string(vy_expand(0,i))
;	;	print, 'Vs,x:'+string(vxs0)
;	;	print, 'Vs,y:'+string(vys0)
;	;	print, 'Vtube,x:'+string(vtubex_should)
;	;	print, 'Vtube,y:'+string(vtubey_should0)
;	;	print, 'Sc 1:'
;	;	print, 'Vobs,x:'+string(vx_expand(1,i))
;	;	print, 'Vobs,y:'+string(vy_expand(1,i))
;	;	print, 'Vs,x:'+string(vxs1)
;	;	print, 'Vs,y:'+string(vys1)
;	;	print, 'Vtube,x:'+string(vtubex_should)
;	;	print, 'Vtube,y:'+string(vtubey_should1)
;	;	stop
;	;endif
;end
;
;i_split = where(abs(dnfront_expand) gt 1.0, n_split)
;if n_split gt 0 then begin
;	print, 'Split events:'
;	print, '# of events: '+string(n_split)
;	print, 'X med: '+string(median(x[*, i_split], /even))
;	print, 'VDFx: '+string(vx_expand[*, i_split])
;	print, 'VDFx med: '+string(median(vx_expand[*, i_split], /even))
;	print, 'Vytube_should: '+string(median(abs(vtubey_should_all[i_split]), /even))
;	print, 'dNy:'
;	print, abs(dnfront_expand(i_split))
;	print, vexpand_should_all(i_split)
;endif else print, 'No split events'
;
;i_head = where((abs(nfront_expand[0,*]) lt 0.2) and (abs(nfront_expand[1,*]) lt 0.2), n_head)
;if n_head gt 0 then begin
;	print, 'Head events:'
;	print, '# of events: '+string(n_head)
;	print, 'X med: '+string(median(x[*, i_head], /even))
;	print, 'Vx med: '+string(median(vx_expand[*, i_head], /even))
;	print, 'Vytube_should all: '+string(abs(vtubey_should_all[i_head]))
;	print, 'Vytube_should: '+string(median(abs(vtubey_should_all[i_head]), /even))
;endif else print, 'No head events'
;
;;; check expansion
;;print, vexpand_should_all
;
;;;; make plots
;;popen, pic_folder+'/dual_dvyex_dny'+t_out_suf
;;print_options,xsize=3.8,ysize=3.
;;stat_plot, abs(dnfront_expand), vexpand_should_all, k_c = 3, bin_range = [0.,2.], binsize = 0.1, qtt_2_range = vy_range, qtt_range = ny_range, qtt_2_title = 'V expand [km/s]', qtt_title = '|'+cap_delta_letter+'n!dy!n'+ny_str+'|', kinbin = kinbin, bin_boundaries = binbounds, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', /no_write, bar_thick = l_thick
;;xyouts, !x.crange(0)+0.5*(!x.crange(1)-!x.crange(0)), !y.crange(1)+0.03*(!y.crange(1)-!y.crange(0)), title, align = 0.5, charsize = 0.9
;;pclose
;;;;;;;;;;;;;;

;;;;;;;;;;;;; various V against ny: fun plot
;vtubey_try = 0. ;in km/s
;
;;vcomp = 'vx'
;;vcomp = 'vy'
;vcomp = 'vtubex'
;
;if strcmp(vcomp, 'vx') then yrange = [min(vx_expand), max(vx_expand)]
;if strcmp(vcomp, 'vy') then yrange = [min(vy_expand), max(vy_expand)]
;if strcmp(vcomp, 'vtubex') then yrange = [-1000., 2000.]
;plot, [-1,1.], yrange, xtitle = 'ny', ytitle = vcomp+' [km/s]', /nodata
;
;for i = 0, n_elements(dnfront_expand)-1 do begin
;	;; probe 0
;	vys0 = vy_expand(0,i)-vtubey_try
;	vxs0 = -vys0*nfront_expand(0,i)/np_expand(0,i)
;	vtube_x0 = vx_expand(0,i)-vxs0
;	;; probe 1
;	vys1 = vy_expand(1,i)-vtubey_try
;	vxs1 = -vys1*nfront_expand(1,i)/np_expand(1,i)
;	vtube_x1 = vx_expand(1,i)-vxs1
;	if strcmp(vcomp, 'vx') then begin
;		oplot, nfront_expand[*, i], vx_expand[*, i]
;		oplot, nfront_expand[*, i], vx_expand[*, i], psym = 4
;	endif
;	if strcmp(vcomp, 'vy') then begin
;		oplot, nfront_expand[*, i], vy_expand[*, i]
;		oplot, nfront_expand[*, i], vx_expand[*, i], psym = 4
;	endif
;	if strcmp(vcomp, 'vtubex') and (n_elements(vtubey_try) lt 2) then begin
;		oplot, nfront_expand[*, i], [vtube_x0, vtube_x1]
;		oplot, nfront_expand[*, i], [vtube_x0, vtube_x1], psym = 4
;	endif
;end
;makepng, pic_folder+'/'+vcomp+'_vs_ny'
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; x dependence of radius and convex/concave ratio ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;; set the range of events
;;xbin_bounds = [-15., -12., -11, -10, -9, -8, -6] ;; highest resolution, used for ratio_convex, 15s.
;xbin_bounds = [-15., -11., -8.5, -6] ;; for radius, 15s 


store_data, 'convex_ratio', data = {xbin_bounds:[-15., -12., -11, -10, -9, -8, -6], x_abc:0.8, y_abc:0.1, yrange:[0.,1.]}
store_data, 'radius', data = {xbin_bounds:[-15., -11., -8.5, -6], x_abc:0.8, y_abc:0.8, yrange:[0., 2.9]}
store_data, 'width', data = {xbin_bounds:[-15., -11., -8.5, -6], x_abc:0.8, y_abc:0.8, yrange:[0., 5.8]}
;things = ['convex_ratio', 'radius']
;things = ['convex_ratio', 'width']
things = ['width']

;;; plot parameters
xrange_plot = [-6., -15.]
abc = ['(b)', '(c)']
left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.013
rate_bar = 0.1
n_panels = n_elements(things)

positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)

popen, pic_folder+'/convex_radius' ;; double plot
;print_options,xsize=2.6,ysize=3.5
print_options,xsize=2.6,ysize=2. ;; single plot
for i_plot = 0, n_panels-1 do begin
	if i_plot eq n_panels-1 then begin
	    xticknames = ''
		xtitle = 'X [R!dE!n]'
	endif else begin
	    xticknames = replicate(' ', 59)
		xtitle = ''
	endelse
	get_data, things(i_plot), data = thing
	xbin_bounds = thing.xbin_bounds
	;; arrays for results
	n_concave_arr = intarr(n_elements(xbin_bounds)-1)
	n_concave_arr(*) = !values.f_nan
	n_convex_arr = n_concave_arr
	n_expand_arr = n_convex_arr
	n_contract_arr = n_convex_arr
	r_mean_arr = fltarr(n_elements(xbin_bounds)-1)+!values.f_nan
	r_std_arr = r_mean_arr
	r_stat = fltarr(n_elements(xbin_bounds)-1, 3)+!values.f_nan ;; lower quartile, median, higher quartile
	dfront_stat = r_stat
	dfront_concave_stat = r_stat
	dfront_convex_stat = r_stat
	dfront_expand_stat = r_stat
	dfront_contract_stat = r_stat
	for i_dis = 0, n_elements(xbin_bounds)-2 do begin
		xrange = xbin_bounds(i_dis:i_dis+1)
	
		i_this = where((x(0,*) gt xrange(0)) and (x(0,*) lt xrange(1)) and (x(1,*) gt xrange(0)) and (x(1,*) lt xrange(1)), j_this)
		if j_this gt 0 then begin
			dfront_this = dfront[i_this]
			dnangle_this = dnangle[i_this]
			rstat, abs(dfront_this), med_dfront, hing1, hing2
			dfront_stat[i_dis, *] = [[hing1], [med_dfront], [hing2]]
		endif
	
		i_this_convex = where((x_convex(0,*) gt xrange(0)) and (x_convex(0,*) lt xrange(1)) and (x_convex(1,*) gt xrange(0)) and (x_convex(1,*) lt xrange(1)), j_this_convex)
		n_convex_arr(i_dis) = j_this_convex
		if j_this_convex gt 0 then begin
			r_this = r_convex(i_this_convex)
			if strcmp(large_dy, 'yes') then begin
				dfront_convex_this = dfront_convex(i_this_convex)
				i_small_dy = where(abs(dfront_convex_this) gt 0.5, n_small_dy)
				if n_small_dy gt 0 then r_bin = r_this(i_small_dy) else r_bin = !values.f_nan
			endif else r_bin = r_this
			rstat, r_bin, med_r, hing1, hing2
;			if i_plot eq 1 then stop
			print, '-- Median radius:'
			print, med_r
			r_stat(i_dis, *) = [[hing1], [med_r], [hing2]]
		endif
	
		i_this_concave = where((x_concave(0,*) gt xrange(0)) and (x_concave(0,*) lt xrange(1)) and (x_concave(1,*) gt xrange(0)) and (x_concave(1,*) lt xrange(1)), j_this_concave)
		n_concave_arr(i_dis) = j_this_concave
	
		i_this_expand = where((x_expand(0,*) gt xrange(0)) and (x_expand(0,*) lt xrange(1)) and (x_expand(1,*) gt xrange(0)) and (x_expand(1,*) lt xrange(1)), j_this_expand)
		n_expand_arr(i_dis) = j_this_expand
	
		i_this_contra = where((x_contra(0,*) gt xrange(0)) and (x_contra(0,*) lt xrange(1)) and (x_contra(1,*) gt xrange(0)) and (x_contra(1,*) lt xrange(1)), j_this_contra)
		n_contract_arr(i_dis) = j_this_contra
	endfor
	;;;;;; plot the variations
	xbincntrs = 0.5*(xbin_bounds(1:*)+xbin_bounds(0:n_elements(xbin_bounds)-2))
	if strcmp(things(i_plot), 'convex_ratio') then begin
		;;; convex ratio
		usersym, 1.*[-1,0,1,0,-1], 1.*[0,1,0,-1,0], thick = l_thick;, /fill ;; diamond
		ratios = double(n_convex_arr)/double(n_convex_arr+n_concave_arr)*100
		plot, xbincntrs, ratios, xrange = xrange_plot, xtitle = xtitle, ytitle = 'Convex/All Percentage', position = positions(*,-(i_plot+1)), /noerase, xtickname = xticknames, xstyle = 1, thick = l_thick-0.4
		oplot, xbincntrs, ratios, psym = 8
		xyouts, xbincntrs, ratios-0.08*(!y.crange(1)-!y.crange(0)), strcompress(string(n_convex_arr), /remove)+'/'+strcompress(string(n_convex_arr+n_concave_arr), /remove), align = 0.5, charsize = 0.5
		;makepng, pic_folder+'/ratio_convex'
	endif
	;;; expand ratio
	;plot, xbincntrs, double(n_expand_arr)/double(n_expand_arr+n_contract_arr), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'Ratio: Expand/Convex'
	;oplot, xbincntrs, double(n_expand_arr)/double(n_expand_arr+n_contract_arr), psym = 4
	;makepng, pic_folder+'/ratio_expand'
	;;;; radius of DFB
	if strcmp(things(i_plot), 'radius') or strcmp(things(i_plot), 'width') then begin
		if strcmp(things(i_plot), 'radius') then begin
			string_q = 'Radius'
			factor = 1.
		endif
		if strcmp(things(i_plot), 'width') then begin
			string_q = 'Width'
			factor = 2.
		endif
		usersym, 0.8*[-1,-1,1,1], 0.8*[1,-1,-1,1], /fill ;; square
		plot, xbincntrs, factor*r_stat(*,1), xrange = xrange_plot, xtitle = xtitle, yrange = thing.yrange, ystyle = 1, ytitle = 'DFB '+string_q+' [R!dE!n]', position = positions(*,n_elements(positions(0,*))-(i_plot+1)), /noerase, xtickname = xticknames, xstyle = 1, /nodata
		oplot, xbincntrs, factor*r_stat(*,1), color = 6, thick = l_thick
		oplot, xbincntrs, factor*r_stat(*,1), psym = 8, color = 6
		oplot, xbincntrs, factor*r_stat(*,0), color = 2
		oplot, xbincntrs, factor*r_stat(*,0), psym = 8, color = 2, symsize = 0.7
		oplot, xbincntrs, factor*r_stat(*,2), color = 2
		oplot, xbincntrs, factor*r_stat(*,2), psym = 8, color = 2, symsize = 0.7
		;makepng, pic_folder+'/r_dfb'
		print, n_convex_arr
		print, total(n_convex_arr)
	endif
	;;; satellite separation in y
	;plot, xbincntrs, dfront_stat(*,1), xrange = [-6., -16.], xtitle = 'X [RE]', ytitle = 'Probes |dY| [RE]', color = 6
	;oplot, xbincntrs, dfront_stat(*,1), psym = 4, color = 6
	;oplot, xbincntrs, dfront_stat(*,0), color = 2
	;oplot, xbincntrs, dfront_stat(*,0), psym = 4, color = 2
	;oplot, xbincntrs, dfront_stat(*,2), color = 2
	;oplot, xbincntrs, dfront_stat(*,2), psym = 4, color = 2
	;makepng, pic_folder+'/dy_probes'
	;;;; plot the bin boundaries
	for i_bound = 0, n_elements(xbin_bounds)-1 do begin 
		if strcmp(things(i_plot), 'radius') or strcmp(things(i_plot), 'width') then y_bars = [max(!y.crange)-1.2*rate_bar*max(!y.crange), max(!y.crange)]
		if strcmp(things(i_plot), 'convex_ratio') then y_bars = [0, 1.2*rate_bar*max(!y.crange)]
		oplot, [xbin_bounds(i_bound), xbin_bounds(i_bound)], y_bars, thick = 0.1
	endfor
	;;;; plot label
	xyouts, !x.crange(0)+thing.x_abc*(!x.crange(1)-!x.crange(0)), !y.crange(0)+thing.y_abc*(!y.crange(1)-!y.crange(0)), abc(i_plot)
endfor
pclose
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; bar plot of dnangle and dy for different type of events to check what controls them ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; settings of binsizes
;binsize_dnangle = 10.
;binsize_dfront = 0.2
;bin_range_angle = [0., 180.]
;bin_range_dfront = [0., 4]
;
;store_data, 'dnangle_convex_concave', data = {good: dnangle_convex, bad: dnangle_concave, bin_range: bin_range_angle, binsize: binsize_dnangle, xtitle: 'dn_angle [degrees]', ytitle: '# of dual-observations', title: 'Convex vs Concave events', xrange: bin_range_angle}
;store_data, 'dnangle_expand_contra', data = {good: dnangle_expand, bad: dnangle_contra, bin_range: bin_range_angle, binsize: binsize_dnangle, xtitle: 'dn_angle [degrees]', ytitle: '# of dual-observations', title: 'Expand vs Contract events', xrange: bin_range_angle}
;store_data, 'dfront_convex_concave', data = {good: dfront_convex, bad: dfront_concave, bin_range: bin_range_dfront, binsize: binsize_dfront, xtitle: 'dfront [RE]', ytitle: '# of dual-observations', title: 'Convex vs Concave events', xrange: bin_range_dfront}
;store_data, 'dfront_expand_contra', data = {good: dfront_expand, bad: dfront_contra, bin_range: bin_range_dfront, binsize: binsize_dfront, xtitle: 'dfront [RE]', ytitle: '# of dual-observations', title: 'Expand vs Contract events', xrange: bin_range_dfront}
;
;groups_check = ['dnangle_convex_concave', 'dnangle_expand_contra', 'dfront_convex_concave', 'dfront_expand_contra']
;
;for i = 0, n_elements(groups_check)-1 do begin 
;	get_data, groups_check(i), data = data
;	bin1d, abs(data.good), abs(data.good), data.bin_range(0), data.bin_range(1), data.binsize(0), kinbin_good, bin_cntrs, flag4nodata = !values.f_nan
;	bin1d, abs(data.bad), abs(data.bad), data.bin_range(0), data.bin_range(1), data.binsize(0), kinbin_bad, bin_cntrs, flag4nodata = !values.f_nan
;	bar_plot_easy, bin_cntrs, kinbin_good+kinbin_bad, xtitle = data.xtitle, ytitle = '# of dual-observations', title=data.title, xrange = data.xrange, xstyle=1
;	bar_plot_easy, bin_cntrs, kinbin_good, /add, color = 6
;	rstat, abs(data.good), good_med, good_lq, good_hq
;	rstat, abs(data.bad), bad_med, bad_lq, bad_hq
;	xyouts, 0.3, 0.7, 'Red: Lowq:'+strcompress(string(good_lq), /remove)+'; Median:'+strcompress(string(good_med), /remove)+'; Highq:'+strcompress(string(good_hq), /remove)+'; Max: '+strcompress(string(max(abs(data.good), /nan)), /remove)+$
;	'!CBlack: Lowq:'+strcompress(string(bad_lq), /remove)+'; Median:'+strcompress(string(bad_med), /remove)+'; Highq:'+strcompress(string(bad_hq), /remove)+'; Max: '+strcompress(string(max(abs(data.bad), /nan)), /remove), /normal
;	makepng, pic_folder+'/'+groups_check(i)+list_suffix
;endfor
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

stop
end
