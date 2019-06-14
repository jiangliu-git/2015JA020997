pro main_estimate_exp
;;; estimate the expansion width and rate of the DFB using data from Liu et al., 2014.
thm_init
@folders

;; data from Liu et al. 2014
x_length = [-20.3, -17.35, -15.2, -13.05, -11.15, -9.75, -8.35, -6.75]
length = [2.68071, 2.01386, 1.86736, 0.975008, 0.789961, 0.519954, 0.402847, 0.307523]
x_bz = [-20.9, -16.5, -14.1, -12.05, -10.3, -8.75, -7.]
bz = [8.31028, 8.27705, 7.98683, 8.61190, 8.91469, 9.03745, 9.54003]
x_vx = [-20.6, -17.625, -15.125, -13., -11.15, -9.75, -8.35, -6.75]
vx = [207.849, 230.204, 306.588, 260.325, 234.247, 180.821, 162.538, 150.550]

x_sep = -12. ;; in RE

width_0 = [0.1, 0.5, 1, 2, 3]
xrange_w = [-22, -6.]
npts = 200
x_width = findgen(npts)/(npts-1)*(xrange_w(1)-xrange_w(0))+xrange_w(0)
dx_width = x_width(1:*)-x_width(0:-2)

length_all = interpol(length, x_length, x_width, /spline)
bz_all = interpol(bz, x_bz, x_width, /spline)
vx_all = interpol(vx, x_vx, x_width, /spline)

;; compute the time difference at two locations for later use of two expansion speeds.
i_near = where(x_width gt x_sep)
i_far = where(x_width lt x_sep)
no_use = min(x_width-x_sep, i_sep, /absolute)
x_near = x_width(i_near)
x_far = x_width(i_far)
vx_near = vx_all(i_near)
vx_far = vx_all(i_far)
dx_near = x_near(1:*)-x_near(0:-2)
dx_far = x_far(1:*)-x_far(0:-2)
dt_near = total(1./vx_near*dx_near)
dt_far = total(1./vx_far*dx_far)

width_mat = dblarr(npts, n_elements(width_0))
vy_exp_mat = dblarr(npts-1, n_elements(width_0))
vy_sep_mat = dblarr(2, n_elements(width_0))

for i = 0, n_elements(width_0)-1 do begin
	flux = width_0(i)*length_all(0)*bz_all(0)	
	width = flux/(length_all*bz_all) ;; plot quantity 1
	dwidth = width(1:*)-width(0:-2)
	vy_exp = vx_all(0:-2)*dwidth/dx_width ;; plot quantity 2
	width_mat[*, i] = width
	vy_exp_mat[*, i] = vy_exp
	vy_near = (width(-1)-width(i_sep))/dt_near
	vy_far = (width(i_sep)-width(0))/dt_near
	vy_sep_mat[*, i] = [vy_near, vy_far]
endfor

;;; make plot
abc = ['(a)', '(b)']
xrange_plot = [xrange_w(1)+1., xrange_w(0)-1.]
;color_arr = [50, 90, 132, 197, 250]
color_arr = [50, 90, 0, 197, 250]

store_data, 'width', data={x:x_width, data:width_mat, range:[0., 20.], title:'DFB Width [R!dE!n]'}
store_data, 'vy', data={x:x_width(0:-2), data:vy_exp_mat, range:[0., 1000.], title:'Exp Rate [km/s]'}
vars = ['width', 'vy']

left_margin = 0.15
right_margin = 0.01
top_margin = 0.05
bot_margin = 0.05
vspace = 0.007
;n_panels = 2
n_panels = 1

positions = panel_positions([1, n_panels], lr_margins = [left_margin, right_margin], bt_margins = [top_margin, bot_margin], space = [0., vspace], height = height)
usersym, 1.*[-1,0,1,0,-1], 1.*[0,1,0,-1,0], /fill ;; diamond

popen, pic_folder+'/estimate'
;print_options,xsize=2.6,ysize=4.3
print_options,xsize=2.6,ysize=2
for i = 0, 0 do begin
;for i = 0, n_panels-1 do begin
	if i eq 0 then title = '' else title = ''
	if i eq n_panels-1 then begin
	    xticknames = ''
		qtt_title = 'X [R!dE!n]'
	endif else begin
	    xticknames = replicate(' ', 59)
		qtt_title = ''
	endelse
	get_data, vars(i), data = this
	plot, this.x, this.data(*,0), /nodata, xrange = xrange_plot, xstyle = 1, yrange = this.range, ystyle = 1, position = positions(*,n_panels-1-i), xtitle = qtt_title, ytitle = this.title, xtickname = xticknames, /noerase, title = title, thick = l_thick
	for j = 0, n_elements(this.data(0,*))-1 do begin
		oplot, this.x, this.data(*,j), color = color_arr(j)
		if j eq 2 then begin
			oplot, this.x[[0, i_sep, n_elements(this.x)-1]], this.data[[0, i_sep, n_elements(this.x)-1], j], color = color_arr(j), line = 1
			oplot, this.x[[0, i_sep, n_elements(this.x)-1]], this.data[[0, i_sep, n_elements(this.x)-1], j], psym = 8, symsize = 0.8, color = color_arr(j)
			print, vy_sep_mat[*,j]
			xyouts, 0.5*(this.x[i_sep]+this.x[-1]), 0.5*(this.data[i_sep,j]+this.data[-1,j])+0.25, '194 km/s', align = 0.5, orientation = -38, /data, charsize = 0.6
			xyouts, 0.5*(this.x[i_sep]+this.x[0]), 0.5*(this.data[i_sep,j]+this.data[0,j])-0.93, '75 km/s', align = 0.5, orientation = -10, /data, charsize = 0.6
		endif
	endfor
	;xyouts, 0.07*(!x.crange[1]-!x.crange[0])+!x.crange[0], 0.85*(!y.crange[1]-!y.crange[0])+!y.crange[0], abc(i)
	if i eq 0 then begin
		xyouts, 0.56, 0.98, 'Estimating DFB expansion', align = 0.5, /normal
	endif
endfor
pclose
stop
end
