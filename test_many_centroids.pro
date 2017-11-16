
pro test_many_centroids,wavelength,n_test,filename=filename,no_cosmic=no_cosmic,nd=nd, $
  method=method,highsky=highsky,req_val=req_val,hsize=hsize,exp_time=exp_time,reads_per_sec=reads_per_sec

;fluxes = [10.0,5.0,3.0,1.0,0.5,0.3,0.1,0.05,0.03,0.01,0.005,0.003,0.001, $
;          0.0005,0.0003,0.00025,0.0002,0.00015,0.00012,0.00011,0.0001,0.00005,0.00003,0.00001]
if (keyword_set(nd)) then begin
    fluxes = [1000000.,500000.,300000.,100000.,50000.,30000.,10000.,5000.,3000.,1000.,500.,300.,100.,50.,30.,10.0,5.0,3.0,1.0,0.5,0.3,0.1]
endif else begin
    fluxes = [10.0,5.0,3.0,1.0,0.5,0.3,0.1,0.05,0.03,0.01,0.005,0.003,0.001, $
              0.0005,0.0003,0.00025,0.0002,0.00015,0.00012,0.00011,0.0001]
endelse
;fluxes = [0.1]
n_fluxes = n_elements(fluxes)

; setup output file
if (not keyword_set(filename)) then filename = 'miri_test_centroids.dat'
openw,unit1,filename,/get_lun
printf,unit1,'# centroids for miri low res faint source tests'

for i = 0,(n_fluxes-1) do begin
    test_centroid,fluxes[i],wavelength,n_test,out_vals=out_vals,/silent,no_cosmic=no_cosmic, $
                  method=method,highsky=highsky,req_val=req_val,hsize=hsize,exp_time=exp_time, $
                  reads_per_sec=reads_per_sec
    printf,unit1,fluxes[i],out_vals[0,*],out_vals[1,*],format='(E10.2,2x,10(E10.2,2x))'
endfor

free_lun,unit1

end
