;dataroot_folder = '/Weiyun/data'
dataroot_folder = '/Weiyun/data/themis'
plasma_folder = '/plasma_new'

if keyword_set(computer) then begin
	fgs_folder=computer + dataroot_folder + '/fgs'
	fgl_folder=computer + dataroot_folder + '/fgl'
	pos_folder=computer + dataroot_folder + '/pos'
	efs_folder=computer + dataroot_folder + '/efs_dsl'
	eff_folder=computer + dataroot_folder + '/eff_dot0_dsl'
	Pall_folder=computer + dataroot_folder + plasma_folder + '/Pall'
	Pth_folder=computer + dataroot_folder + plasma_folder + '/Pth'
	vi_folder = computer + dataroot_folder + plasma_folder + '/vi'
	viperp_folder = computer + dataroot_folder + plasma_folder + '/viperp'
	ve_folder = computer + dataroot_folder + plasma_folder + '/ve'
	vexb_folder = computer + dataroot_folder + plasma_folder + '/vexb'
	veperp_folder = computer + dataroot_folder + plasma_folder + '/veperp'
	Pttl_fgs_folder=computer+dataroot_folder+plasma_folder+'/Pttl_fgs'
	Blobe_folder=computer+dataroot_folder+plasma_folder+'/Blobe'
endif

listfolder='../lists'
save_folder='variables'
pic_folder = '../../../temp'

;;;;;; 
l_thick = 2.4
;;; constants and signs
Re = 6371.
mu0 = 4*!pi*1e-7
nTesla2_to_nPa = 0.01/25.132741
perp_sign = '!9'+string("136B)+'!X'
cross = '!9'+string("264B)+'!X'
theta_letter = '!9'+string("161B)+'!X'
cap_theta_letter = '!9'+string("121B)+'!X'
delta_letter = '!9'+string("144B)+'!X'
cap_delta_letter = '!9'+string("104B)+'!X'

;;; for DFB properties
seconds_check = 15.

;;; turn off the timestamp of tplot
time_stamp,/off
