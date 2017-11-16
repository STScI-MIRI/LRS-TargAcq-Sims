
pro fit_exptime_flux_miri,wavelength,req_val,prompt=prompt,eps=eps,subarray_exp=subarray_exp

exptime_str = ['12','120','240','360','480','720','960']
files = 'miri_run_h2_' + wavelength + '_m2_r1_exptime' + $
        exptime_str + '.dat'
ps_file = 'miri_run_h2_' + wavelength + '_m2_r1_exptime_fit'
n_files = n_elements(files)

exptimes = float(exptime_str)
est_req_fluxes = fltarr(n_files)
for i = 0,(n_files-1) do begin
    plot_many_centroids,files[i],/radial,req_val=req_val,est_req_flux=est_req_flux,eps=eps
    est_req_fluxes[i] = est_req_flux
    if (keyword_set(prompt)) then begin
        tstr = ''
        read,'Continue: ',tstr
    endif
endfor

if (keyword_set(subarray_exp)) then begin
    est_req_fluxes *= sqrt(30./subarray_exp)
    exptimes *= subarray_exp/30.
endif

loadct,0
fit_exptime_flux,exptimes,est_req_fluxes,eps=eps,ps_file=ps_file,req_val=req_val

print,est_req_fluxes

; F560W, MRS-IFU, 20 mas
;fit_exptime_flux,[120.,240.,360.,480.,720.,960.],[4.,2.6,1.6,0.9,0.7,0.6]

end

