; program to create a MIRI image w/ the ramps
;  input source flux in mJy, wavelength in microns

pro create_miri_image,source_flux,wavelength,show_plot=show_plot,seed=seed, $
  image=image,est_image=est_image,sky_bkg_out=sky_bkg_out,exp_center=exp_center, $
  write_images=write_images,shift_x=shift_x,shift_y=shift_y,silent=silent,no_cosmic=no_cosmic, $
  method=method,highsky=highsky,exp_time=exp_time,reads_per_sec=reads_per_sec,subarray=subarray

if (n_params() LT 1) then begin
    print,'create_miri_image,source_flux,wavelength'
    print,'  source_flux in mJy and wavelength in microns'
    return
endif

if (not keyword_set(method)) then method = 1

in_source_flux = source_flux

n_rej = 1
;total_eff = 0.25
if (not keyword_set(exp_time)) then exp_time = 960.
if (not keyword_set(reads_per_sec)) then reads_per_sec = 1./30.
; make sure the # of reads is 4n+1
tn_reads = round(exp_time*reads_per_sec + 1)
tn_sub = round((tn_reads - 1)/4.)
nexp_time = (tn_sub*4.)/reads_per_sec
if (nexp_time NE exp_time) then begin
    exp_time = nexp_time
    print,'new exptime (4n+1) = ', exp_time
endif

gain = 5.96
read_noise = 19. ; per 30s read
if (keyword_set(subarray)) then begin
    if (reads_per_sec GT 1./1.) then read_noise *= sqrt(8.)
endif else begin
    if (reads_per_sec GT 1./29.) then read_noise *= sqrt(8.)
endelse
;read_noise = read_noise/sqrt(8) ; per 30s read (average of 8 reads)
;read_noise /= sqrt(12./34.) ; correcting to a single 30 sec read
dark_current = 0.3 ; e/sec
;sky_bkg = 50. ; photons/sec/micron/pixel
;sky_bkg *= total_eff  ; accounting for system throughput/QE

cr_prob_per_sec_per_pixel = 0.5/1000.       ;  50% of pixels affected in 1000s
cr_prob_per_int_per_pixel = cr_prob_per_sec_per_pixel*exp_time

; convert the input source flux from mJy to photons/sec
;source_flux *= 1h-26*1e-3  ; from mJy to W m^-2 Hz^-1
;source_flux /= 6.626d-34*(3.d14/wavelength)  ; convert to photons m^-2 Hz^-1
;source_flux *= ; multiply by bandwidth in Hz

; conversion factor for G. Rieke's spreadsheet (depends on wavelength)
; updated 4 Apr 2008 to reflect latest spreadsheet MIRIradmodeldb.xls
; (same as da.xls)
case wavelength of
    5.6 : begin
        conv_factor = 1.26e4
        sky_bkg = 7.46  ; in e/sec
        fwhm = 0.2
    end
    10.0 : begin
        conv_factor = 1.62e4
        sky_bkg = 93.7  ; in e/sec
        fwhm = 0.35
    end
    13.5 : begin
        conv_factor = 27.1
        sky_bkg = 3.21  ; in e/sec
        fwhm = 0.47
    end
    15.0 : begin
        conv_factor = 1.34e4
        sky_bkg = 366.  ; in e/sec
        fwhm = 0.53
    end
    25.5 : begin
        conv_factor = 8.5e3
        sky_bkg = 8280.  ; in e/sec
        fwhm = 0.9
    end
    else : stop
endcase
source_flux *= conv_factor/gain
pixscale = 0.11

if (keyword_set(highsky)) then begin
    sky_bkg *= 2.
endif
sky_bkg_out = (sky_bkg/reads_per_sec)/gain

image = get_miri_pixelated_image(fwhm/pixscale,shift_x,shift_y,100,32,exp_center,write_images=write_images)
image_size = size(image)

; determine the source e/sec  ; make it 1 e/s max
total_psf = total(image)
image *= (source_flux/total_psf)
;print,'max source flux in central pixel = ', max(image)

; add the background (including dark current)
image += sky_bkg/gain + dark_current/gain

; write out noiseless image
if (keyword_set(write_images)) then begin
    fits_write,'miri_no_noise.fits',image,header
endif

; crop the image
;hsize = 16
;x1 = image_size[1]/2 - hsize
;x2 = image_size[1]/2 + hsize
;y1 = image_size[2]/2 - hsize
;y2 = image_size[2]/2 + hsize
;image = image[x1:x2,y1:y2]
;image_size = size(image)

; generate a zero_point image
zero_point = replicate(5000.0,image_size[1],image_size[2])

; loop over the image creating ramps for each pixel
;   image in e/sec
fit_slope = fltarr(image_size[1],image_size[2])
fit_slope2 = fltarr(image_size[1],image_size[2])
fit_slope_unc = fltarr(image_size[1],image_size[2])
fit_zeropt = fltarr(image_size[1],image_size[2])
n_cr = 0
for i = 0,(image_size[1]-1) do begin
;    print,i
    for j = 0,(image_size[2]-1) do begin
        create_ramp,zero_point[i,j],image[i,j],exp_time,reads,dn_vals,dn_unc_vals, $
                    dn_ran_unc_vals,dn_cor_unc_vals, $
                    read_noise=read_noise,gain=gain,reads_per_sec=reads_per_sec, $
                    seed=seed

        ; remove the first read (which is at
        ; the start of the ramp - no time elapsed
        n_dn_vals = n_elements(dn_vals)
        reads = reads[1:n_dn_vals-1]
        dn_vals = dn_vals[1:n_dn_vals-1]
        dn_unc_vals = dn_unc_vals[1:n_dn_vals-1]

        n_dn_vals = n_elements(dn_vals)
        if (keyword_set(show_plot)) then begin
            print,dn_vals
            kplot,reads,dn_vals,psym=1
            ans = ''
            read,'Continue: ',ans
        endif

        ; add cosmic rays if asked
        if (not keyword_set(no_cosmic)) then begin
            val = randomu(seed)
            if (val LE cr_prob_per_int_per_pixel) then begin
                read_num = fix(randomu(seed)*n_dn_vals)
                cr_amp = 20000./gain   ; in DN
                dn_vals[read_num:n_dn_vals-1] += cr_amp
                n_cr++
            endif
        endif

        if ((i+j) EQ 0) then begin
            n_dn_vals = n_elements(dn_vals)
            cube = fltarr(image_size[1],image_size[2],n_dn_vals)
        endif
        cube[i,j,*] = dn_vals

        ; find the cosmic ray (if it exists)
        diffs = dn_vals[n_rej+1:n_dn_vals-1] - dn_vals[n_rej:n_dn_vals-2]
        indxs = indgen(n_elements(diffs))
        if (n_elements(diffs) GT 4) then begin
            iter_sigma_clip,diffs,indxs,sigma_factor=5.,sigma_vals=sigma_vals,silent=1
        endif
        if (n_elements(indxs) LT n_elements(diffs)) then begin
            indxs = where(abs(diffs-sigma_vals[0]) GT 5.*sigma_vals[2],n_indxs)
            n_cr = min(indxs) + n_rej
            if (n_cr GT n_dn_vals/2.) then begin
                n_beg = n_rej
                n_end = n_cr
            endif else begin
                n_beg = n_cr + 1
                n_end = n_dn_vals - 1
            endelse
        endif else begin
            n_beg = n_rej
            n_end = n_dn_vals-1
        endelse
        ; fit ramps and derive the slope image
        if (n_end - n_rej GE 2) then begin
            fit_coeff = kg_linfit_cor(reads[n_beg:n_end],dn_vals[n_beg:n_end], $
                                      dn_unc_vals[n_beg:n_end],read_noise, $
                                      yfit=yfit,sigma=fit_coeff_unc)
        endif else begin
;            print,'***crap***'
            fit_coeff = [0.0,0.0]
            fit_coeff_unc = [0.0,0.0]
;            stop
        endelse
        fit_slope[i,j] = fit_coeff[1]
        fit_slope_unc[i,j] = fit_coeff_unc[1]
        fit_zeropt[i,j] = fit_coeff[0]

        ; create the slope image in a different way
        ; this is the current way the flight software would do it
        
        ; first average down to 3 images (coadd)
        if (method EQ 1) then begin
            dn_vals_rebin = rebin(dn_vals,4)
            reads_rebin = rebin(reads,4)
            slope_vals = ([dn_vals_rebin[2:3]] - [dn_vals_rebin[1:2]])/((n_dn_vals/4.)*[1.,1.]);
            slope_vals = [1.0]
        endif else if (method EQ 2) then begin
            ; make 2 2pt differences (updated slope estimation)
            bindxs = (indgen(4))*fix((n_dn_vals/4.))
            eindxs = (indgen(4) + 1)*fix((n_dn_vals/4.)) - 1
            slope_vals = (dn_vals[eindxs] - dn_vals[bindxs])/(n_dn_vals/4.)
        endif else begin
            if (n_elements(reads) NE 4) then begin
                print,'argh -> n_reads = ', n_elements(reads)
                stop
            endif
            slope_vals = ([dn_vals[2:3]] - [dn_vals[1:2]])/((n_dn_vals/4.)*[1.,1.])
        endelse
        fit_slope2[i,j] = min(slope_vals)
    endfor
endfor

;stop

if (keyword_set(write_images)) then begin
; write image
    fits_write,'miri_test.fits',cube,header

; write reduced images
    fits_write,'miri_slope.fits',fit_slope,header
    fits_write,'miri_slope_unc.fits',fit_slope_unc,header
    fits_write,'miri_zeropt.fits',fit_zeropt,header

    fits_write,'miri_slope2.fits',fit_slope2,header
endif

if (not keyword_set(silent)) then print,'n_cr = ',n_cr

source_flux = in_source_flux
image = fit_slope
est_image = fit_slope2

end
