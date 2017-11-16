; program to create a bunch of MIRI images and see how well the
; centroid can be determined
; written by KDG 11 Oct 2007 to support the testing for the
; faint peakup for the low resolution spectrometer

pro test_centroid,source_flux,wavelength,n_test,write_images=write_images,method=method, $
                  out_vals=out_vals,use_input_sky=use_input_sky,silent=silent,no_cosmic=no_cosmic, $
                  highsky=highsky,req_val=req_val,hsize=hsize,exp_time=exp_time,reads_per_sec=reads_per_sec

iter_thres = 0.01
max_iter = 15
if (not keyword_set(req_val)) then req_val = 1.
if (not keyword_set(hsize)) then hsize = 2

;cor_center = [16.,16.]
case wavelength of
    5.6 : begin
        zero_off = 0.0078 + (-1*[-0.495,-0.495]) + [-2.00644e-05,0.000449146]
    end
    10.0 : begin
        zero_off = [0.496,0.490] + [-0.001666,0.002386]
    end
    13.5 : begin
        zero_off = [0.4755,0.4727] + [-8.17E-02,-1.30E-01]
        hsize = 2
    end
    15.0 : begin
        zero_off = [0.4755,0.4727] + [-8.17E-02,-1.30E-01]
        hsize = 2
    end
    25.5 : begin
        zero_off = [0.849,0.847] + [0.0056,0.00036]
        hsize = 4
    end
    else : stop
endcase
;zero_off = [0.,0.]

save_m1 = fltarr(2,n_test)
save_em1 = fltarr(2,n_test)
for i = 0,(n_test-1) do begin
    ; determine the subpixel offsets
    shift_x = (randomu(seed)*2.) - 1.
    shift_y = (randomu(seed)*2.) - 1.
;    shift_x = 0.
;    shift_y = 0.
    if (not keyword_set(silent)) then begin
        print,'shift values (pixels, x,y) = ', shift_x,shift_y
    endif
    ; make the image
    create_miri_image,source_flux,wavelength,seed=seed,show_plot=0,method=method, $
                      image=image,est_image=est_image,silent=silent,highsky=highsky, $
                      write_images=write_images,sky_bkg_out=sky_bkg_out, $
                      shift_x=shift_x,shift_y=shift_y,exp_center=exp_center,no_cosmic=no_cosmic, $
                      exp_time=exp_time,reads_per_sec=reads_per_sec
    image_size = size(image)

    ; get background
    bkg = median(image)
    bkg2 = mean(image)
    sindxs = sort(image)
    n_sindxs = n_elements(sindxs)
    bkg3 = image[sindxs[fix(n_sindxs*0.3)]]
    if (keyword_set(use_input_sky)) then begin
        bkg = sky_bkg_out
    endif
;    print,'bkg = ', bkg, bkg2, bkg3
    image -= bkg3

    ; get background
    est_bkg = median(est_image)
    est_bkg2 = mean(est_image)
    sindxs = sort(est_image)
    n_sindxs = n_elements(sindxs)
    est_bkg3 = est_image[sindxs[fix(n_sindxs*0.3)]]
    if (keyword_set(use_input_sky)) then begin
        est_bkg = sky_bkg_out
    endif
    est_image -= est_bkg3

    ; crop the image
;    x1 = image_size[1]/2 - hsize
;    x2 = image_size[1]/2 + hsize
;    y1 = image_size[2]/2 - hsize
;    y2 = image_size[2]/2 + hsize
;    sub_image = image[x1:x2,y1:y2] - bkg
;    sub_image_est = est_image[x1:x2,y1:y2] - bkg_est
    find_peak,image,sub_image,hsize,peak_x,peak_y
;    peak_x = fix(round(exp_center[0]))
;    peak_y = fix(round(exp_center[1]))
    if (not keyword_set(silent)) then begin
        print,'rough (find_peak) center = ',peak_x,peak_y
    endif
;    print,sub_image

    find_varwindow_moment,sub_image,iter_thres,max_iter,m1,m2,m3,skew
    save_m1[*,i] = m1 + [peak_x,peak_y] - hsize - exp_center - zero_off;cor_center - [shift_x,shift_y] - zero_off

    find_peak,est_image,sub_image_est,hsize,epeak_x,epeak_y
    if (not keyword_set(silent)) then begin
        print,'rough (efind_peak) center = ',peak_x,peak_y
    endif

    find_varwindow_moment,sub_image_est,iter_thres,max_iter,em1,em2,em3,eskew

    if (not keyword_set(silent)) then begin
        print,'expected center = ',exp_center
    endif
    save_em1[*,i] = em1 + [epeak_x,epeak_y] - hsize - exp_center - zero_off;cor_center - [shift_x,shift_y] - zero_off
    if (not keyword_set(silent)) then print,i,save_m1[*,i],save_em1[*,i]        ;,m2,m3,skew,iter_thres
;    print,i,em1;,m2,m3,skew,iter_thres
;    ans = ''
;    read,ans
endfor

indxs = where(finite(save_m1[0,*]) AND finite(save_m1[1,*]),n_indxs)
if (n_indxs GT 3) then begin
;    image_statistics,save_m1[0,indxs],mean=mean_m1x,stddev=stddev_m1x
;    image_statistics,save_m1[1,indxs],mean=mean_m1y,stddev=stddev_m1y
    mean_m1x = 0.
    mean_m1y = 0.
    indxs = where(finite(save_m1[0,*]) AND finite(save_m1[1,*]) AND $
                  (abs(save_m1[0,*] - mean_m1x) LT req_val) AND (abs(save_m1[1,*] - mean_m1y) LT req_val),n_indxs)
endif
if (n_indxs GT 5) then begin
    image_statistics,save_m1[0,indxs],mean=mean_m1x,stddev=stddev_m1x
    image_statistics,save_m1[1,indxs],mean=mean_m1y,stddev=stddev_m1y

    iter_sigma_clip,save_m1[0,*],indxs,sigma_vals=sigma_vals
    mean_m1x = sigma_vals[0]
    stddev_m1x = sigma_vals[2]
    iter_sigma_clip,save_m1[1,*],indxs,sigma_vals=sigma_vals
    mean_m1y = sigma_vals[0]
    stddev_m1y = sigma_vals[2]
endif else begin
    mean_m1x = 0.
    stddev_m1x = 0.
    mean_m1y = 0.
    stddev_m1y = 0.
endelse

indxs = where(finite(save_em1[0,*]) AND finite(save_em1[1,*]),n_eindxs)
if (n_eindxs GT 3) then begin
;    image_statistics,save_em1[0,indxs],mean=mean_em1x,stddev=stddev_em1x
;    image_statistics,save_em1[1,indxs],mean=mean_em1y,stddev=stddev_em1y
    mean_em1x = 0.
    mean_em1y = 0.
    indxs = where(finite(save_em1[0,*]) AND finite(save_em1[1,*]) AND $
                  (abs(save_em1[0,*] - mean_em1x) LT req_val) AND (abs(save_em1[1,*] - mean_em1y) LT req_val),n_eindxs)
endif
if (n_eindxs GT 5) then begin
    image_statistics,save_em1[0,indxs],mean=mean_em1x,stddev=stddev_em1x
    image_statistics,save_em1[1,indxs],mean=mean_em1y,stddev=stddev_em1y

    iter_sigma_clip,save_em1[0,*],indxs,sigma_vals=sigma_vals
    mean_em1x = sigma_vals[0]
    stddev_em1x = sigma_vals[2]
    iter_sigma_clip,save_em1[1,*],indxs,sigma_vals=sigma_vals
    mean_em1y = sigma_vals[0]
    stddev_em1y = sigma_vals[2]
endif else begin
    mean_em1x = 0.
    stddev_em1x = 0.
    mean_em1y = 0.
    stddev_em1y = 0.
endelse

print,'x = ', mean_m1x, stddev_m1x,n_indxs
print,'ex = ', mean_em1x, stddev_em1x,n_eindxs

print,'y = ', mean_m1y, stddev_m1y
print,'ey = ', mean_em1y, stddev_em1y

out_vals = fltarr(2,5)
out_vals[0,*] = [mean_m1x,stddev_m1x,mean_m1y,stddev_m1y,n_indxs]
out_vals[1,*] = [mean_em1x,stddev_em1x,mean_em1y,stddev_em1y,n_eindxs]

end
