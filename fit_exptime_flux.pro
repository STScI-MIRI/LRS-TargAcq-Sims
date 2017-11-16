
pro fit_exptime_flux,exptime,flux,eps=eps,ps_file=ps_file,req_val=req_val

if (not keyword_set(ps_file)) then ps_file = 'fit_exptime_flux'
setup_ps_output,ps_file,eps=eps

kplot,flux,exptime,psym=1,kplot_type='oo',xtitle='flux [!4l!3Jy]',ytitle='exptime [sec]'

fit_coeff = poly_fit(alog10(flux),alog10(exptime),1,yfit=yfit)

koplot,flux,10^yfit,psym=100,linestyle=2

fit_eq = 'exptime = ' + strtrim(string(fit_coeff[0],format='(F8.2)'),2) + $
         '(flux)!U' + strtrim(string(fit_coeff[1],format='(F8.2)'),2) + '!D'

xyouts,0.5,0.8,fit_eq,/normal,charsize=1.5
if (keyword_set(req_val)) then begin
    xyouts,0.5,0.75,'req. precision = ' + strtrim(string(req_val,format='(F8.1)'),2) + ' mas', $
           /normal,charsize=1.5
endif

print,transpose(fit_coeff)

close_ps_output,eps=eps

end
