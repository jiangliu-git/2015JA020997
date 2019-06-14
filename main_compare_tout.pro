pro main_compare_tout
@folders

method_suffix = 'binxbout'
list_suffix = '_earthward_df_'+method_suffix

;;; load data
t_out_0 = datain_simple(save_folder+'/tout_dual'+list_suffix+'_0'+'.dat', dim=1, type='double')
t_0_0 = datain_simple(save_folder+'/t0_dual'+list_suffix+'_0'+'.dat', dim=1, type='double')
t_out_15 = datain_simple(save_folder+'/tout_dual'+list_suffix+'.dat', dim=1, type='double')
t_0_15 = datain_simple(save_folder+'/t0_dual'+list_suffix+'.dat', dim=1, type='double')
t_out_60 = datain_simple(save_folder+'/tout_dual'+list_suffix+'_60'+'.dat', dim=1, type='double')

;pm, t_0_15-t_out_15
pm, t_0_0-t_out_0

stop
end
