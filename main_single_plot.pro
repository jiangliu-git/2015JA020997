pro main_single_plot

thm_init

computer = 'I:'
;computer = '/home/jliu'

@folders

method_suffix = 'binxbout'

;;; electric field used for vdHT and vperp2
etype = 'eff' 
;etype = 'efs'

;;; time range to get t_out (minus this seconds)
;t_out_suf = '' ;; default: 15s
t_out_suf = '_0'
;t_out_suf = '_3'
;t_out_suf = '_5'

list_suffix = '_earthward_df_'+method_suffix+t_out_suf

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; bin plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; for near and far
x_sep = -12.
y_dusk = 2.
y_dawn = 0.

k_c = 5
bin_range = [-1., 1.]
binsize = 0.2

ny_near_bins = [0.13, 0.39, 0.65, 1.]
ny_far_bins = [0.2, 0.6, 1.]

ny_range = [-1., 1.]
vy_range = [-250., 250.]
vx_range = [-100., 300.]

;;; load data
pos = datain_simple(save_folder+'/pos'+list_suffix+'.dat', dim=3, type='double')
n = datain_simple(save_folder+'/n'+list_suffix+'.dat', dim=3, type='double')
vc = datain_simple(save_folder+'/vperp2_'+etype+list_suffix+'.dat', dim=3, type='float')

;;; components
x = pos(0,*)/RE
y = pos(1,*)/RE
z = pos(2,*)/RE
nx = n[0,*]
ny = n[1,*]
nz = n[2,*]
vc_x = vc(0,*)
vc_y = vc(1,*)
vc_z = vc(2,*)
nxy = n[0:1,*]
nxy_length = transpose(sqrt(total(nxy^2, 1)))
nxy = nxy/[nxy_length, nxy_length]
ny_star = nxy[1,*]

;;; select the type of ny to use
ny_use = ny & ny_title = 'n!dy' 
;ny_use = ny_star & ny_title = 'n!s!dy!r!u*!n'

;;; locations
near_x = x gt x_sep
far_x = x lt x_sep
dusk = y gt y_dusk
dawn = y lt y_dawn

store_data, 'gp_near_dusk', data = {i_good: where(near_x and dusk)}
store_data, 'gp_near_dawn', data = {i_good: where(near_x and dawn)}
store_data, 'gp_near',  data = {i_good: where(near_x), bin_bounds_ny: [-reverse(ny_near_bins), ny_near_bins], loc_str: 'X>-12 R!dE!n'}
store_data, 'gp_far',  data = {i_good: where(far_x), bin_bounds_ny: [-reverse(ny_far_bins), ny_far_bins], loc_str: 'X<-12 R!dE!n'}

;groups = ['gp_near_dusk', 'gp_near_dawn', 'gp_far']
groups = ['gp_near', 'gp_far']

;;; plot parameters
if strcmp(t_out_suf, '_0') then abc = ['(a)', '(b)', '(c)'] else abc = ['(h)', '(i)', '(j)']
left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.007
n_panels = n_elements(groups)

positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)

popen, pic_folder+'/vperp2_bin_'+etype+list_suffix
print_options,xsize=2.6,ysize=4.3
for i = 0, n_panels-1 do begin
	if i eq 0 then title = '' else title = ''
	if i eq n_panels-1 then begin
	    xticknames = ''
		qtt_title = ny_title
	endif else begin
	    xticknames = replicate(' ', 59)
		qtt_title = ''
	endelse
	get_data, groups(i), data = gp
	;stat_plot, ny_use(gp.i_good), vc_y(gp.i_good), k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = vy_range, qtt_range = ny_range, qtt_2_title = 'V!dDF,y!n [km/s]', qtt_title = qtt_title, kinbin = kinbin, bin_boundaries = gp.bin_bounds_ny, title = title, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', position = positions(*,-(i+1)), /noerase, qtt_tickname = xticknames, /no_write
	stat_plot, ny_use(gp.i_good), vc_x(gp.i_good), k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = vx_range, qtt_range = ny_range, qtt_2_title = 'V!dDF,x!n [km/s]', qtt_title = qtt_title, kinbin = kinbin, bin_boundaries = gp.bin_bounds_ny, title = title, /no_mean, color_med = 6, color_quar = 2, type_med = 'square', position = positions(*,-(i+1)), /noerase, qtt_tickname = xticknames, /no_write
	oplot, !x.crange, [0,0]
	oplot, [0,0], !y.crange
	xyouts, 0.07*(!x.crange[1]-!x.crange[0])+!x.crange[0], 0.85*(!y.crange[1]-!y.crange[0])+!y.crange[0], abc(i)
	xyouts, 0.6*(!x.crange[1]-!x.crange[0])+!x.crange[0], 0.06*(!y.crange[1]-!y.crange[0])+!y.crange[0], gp.loc_str
endfor
pclose

;;;;;; general use: modify for binned values 
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
;bin_boundaries = [-24, -17.8, -15.2,-13,-11.1,-9.5, -8,-6]
;
;;;;;; qtt 2
;qtt_2_bin = vc_y
;qtt_2_title = 'Vy [km/s]'
;qtt_2_str = 'vy'
;qtt_2_range = [-250., 250.]
;k_c = 5
;pm = 1
;
;;qtt_2_bin = abs(vc_y)
;;qtt_2_title = '|Vy| [km/s]'
;;qtt_2_str = 'vyabs'
;;qtt_2_range = [0., 250.]
;;k_c = 5
;
;stat_plot, transpose(qtt_bin), transpose(qtt_2_bin), k_c = k_c, bin_range = bin_range, binsize = binsize, qtt_2_range = qtt_2_range, qtt_range = qtt_range, qtt_2_title = qtt_2_title, qtt_title = qtt_title, kinbin = kinbin, bincntrs_out = bincenters, pm = pm, ratio_pm = ratio_pm, vertical = vertical, bin_boundaries = bin_boundaries, title = title, avrg = avrg, std = std, med = med, /no_mean, color_med = 6, color_quar = 2, type_med = 'square'
;oplot, !x.crange, [0,0]
;oplot, [0,0], !y.crange
;makepng, pic_folder+'/'+qtt_2_str+'_'+qtt_str
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;; superpose plot ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;; load tplot variables from saved data and add dashed lines
;
;btype = 'fgl'
;etype = 'efs'
;
;var_names = ['bz_near_'+btype, 'vy_near_evening_'+etype, 'vy_near_morning_'+etype, 'bz_far_'+btype, 'vy_far_evening_'+etype, 'vy_far_morning_'+etype]+list_suffix
;
;for i = 0, n_elements(var_names)-1 do begin
;  name = var_names(i)
;  datain, name, filename = save_folder+'/superpose_'+name+'.dat', dim = 3
;  ;;; need to split vec, inorder to make dash-solid-dash
;  split_vec, name
;  del_data, name
;  options, name+'_x', linestyle = 1
;  options, name+'_z', linestyle = 1
;  store_data, name, data = [name+'_x',name+'_y',name+'_z']
;endfor
;
;options, '*', xtickname=['-1','0','1'], thick=l_thick
;options, 'bz**', ytitle = 'B!dz!n!c', ysubtitle = '[nT]'
;options, 'vy*evening*', ytitle = 'V!dExB,y!n Evening!c', ysubtitle = '[km/s]'
;options, 'vy*morning*', ytitle = 'V!dExB,y!n Morning!c', ysubtitle = '[km/s]'
;;ylim, 'vy*evening', -49., 181.
;;ylim, 'vy*morning', -119., 49.
;
;tplot_options,'vtitle',''
;;;;;;;;;; vy plot
;popen, pic_folder+'/superpose'+list_suffix
;print_options,xsize=5.6,ysize=8.5
;tplot, var_names, $
;trange = ['1999 12 31 23 59', '2000 1 1 0 1'], title = ''
;;; make lines in left and right to indicate what quantity is being plotted.
;timebar, '2000 1 1 0 0'
;for i=0, n_elements(var_names)-1 do begin
;  this_var = var_names(i)
;  timebar, 0, /databar, varname=this_var, line = 3
;endfor
;xyouts, 0.43, 0.03, 'Minutes to t!d0!n', /normal
;
;;;; right margin signs
;hpos = 0.85
;vert = 0.38
;hori = 0.03
;;; Near sign
;pos = [hpos, 0.75]
;make_px, pos, ver=vert, hor=hori, charname='X>-12 R!dE!n', lchar=0.15, orientation = -90, alignment = 0.5
;;; Far sign
;pos = [hpos, 0.3]
;make_px, pos, ver=vert, hor=hori, charname='X<-12 R!dE!n', lchar=0.15, orientation = -90, alignment = 0.5
;
;;;; label panels
;xs=0.21*(1+dblarr(6,1))
;ys=[0.14, 0.34,0.45,0.61,0.76,0.938]
;ss=['(g)', '(f)','(e)','(d)','(c)','(b)']
;xyouts, xs, ys, ss, charsize=1.2, /normal
;;xyouts, xs+0.05, ys, [n_dawn_binxbout, n_dusk_binxbout, n_dawn_binxl, n_dusk_binxl, n_dawn_minvar, n_dusk_minvar], charsize=1, /normal
;pclose

end
