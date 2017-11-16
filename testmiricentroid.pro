pro testmiricentroid,num_subpixel,total_counts,read_noise,expname,num_repeats,poisson_noise=poisson_noise,chisqr_thres=chisqr_thres, $
                 second_moment_thres=second_moment_thres, third_moment_thres=third_moment_thres,bad_pixel_prob=bad_pixel_prob,bad_pixel_value=bad_pixel_value,background=background,use_pixel_mask=use_pixel_mask,win_size=win_size,$
                 random_psf=random_psf,num_links=num_links,fwhm=fwhm,iter_thres=iter_thres,threshold=threshold

secondmoment_threshold = 0.9
thirdmoment_threshold = 0.7
pixel_size = 25e-6              ; pixel size in meters
unshifted_peak=[513,511]
seed=2345678
scaling_factor = 0.0

openw, lun_out, expname+'_errorinfo.txt', /get_lun


if (not keyword_set(threshold)) then begin 
    threshold=-100
endif else begin
    threshold=threshold
endelse

if (not keyword_set(fwhm)) then begin 
    fwhm = 2.8
endif else begin
    fwhm=fwhm
endelse

if (not keyword_set(iter_thres)) then begin 
    iter_thres=0.01
endif else begin
    iter_thres=iter_thres
endelse

if (not keyword_set(chisqr_thres)) then begin 
    chi_squared_threshold = 3.0
endif else begin
    chi_squared_threshold = chisqr_thres
endelse

if (not keyword_set(win_size))then begin
    window_size=2
endif else begin
    window_size=win_size
endelse

if (not keyword_set(second_moment_thres)) then begin 
    secondmoment_threshold = 0.5
endif else begin
    secondmoment_threshold = second_moment_thres
endelse

if (not keyword_set(third_moment_thres)) then begin 
    third_threshold = 0.3
endif else begin
    thirdmoment_threshold = third_moment_thres
endelse

if (not keyword_set(bad_pixel_prob)) then begin
    bad_pixel_probability = 0.0
endif else begin
    bad_pixel_probability = bad_pixel_prob
endelse

if (not keyword_set(bad_pixel_value)) then begin
    bad_value = 0.0
endif else begin
    bad_value = bad_pixel_value
endelse

if (not keyword_set(background)) then begin
    background = 0.0
endif 
if (not keyword_set(use_pixel_mask)) then begin
    use_pixel_mask = 0
endif 

if (not keyword_set(num_links)) then begin
    num_links=0
endif else begin
    num_links=num_links
endelse

model_scale = 1.25e-6
scale_factor = round(pixel_size/model_scale)
binned_peak = unshifted_peak/scale_factor

stepsize = scale_factor/num_subpixel ;
; Either create or read the model cube
xrms = fltarr(11,22)
yrms = fltarr(11,22)
xavg = fltarr(11,22)
yavg = fltarr(11,22)
xmoment2 = fltarr(11,22)
ymoment2 = fltarr(11,22)
xmoment3 = fltarr(11,22)
ymoment3 = fltarr(11,22)
residavg = fltarr(11,22)
xrmsfit  = fltarr(11,22)
yrmsfit  = fltarr(11,22)
xrms_thres  = fltarr(11,22)
yrms_thres  = fltarr(11,22)
fitxrms_thres  = fltarr(11,22)
fityrms_thres  = fltarr(11,22)
centhres_count = intarr(11,22)
fitthres_count = intarr(11,22)
        if (not keyword_set(random_psf)) then begin
            filename='miri_ND_psf.fits'
            fits_read,filename,inimage
        endif
        xsum=0
        x2sum=0
        x2sum_cor=0
        ysum =0
        y2sum =0
        y2sum_cor =0
        varxsum=0
        varx2sum=0
        varysum =0
        vary2sum =0
        gwxsum=0
        gwx2sum=0
        gwysum =0
        gwy2sum =0
        count =0
        moment2sum =fltarr(2)
        moment3sum =fltarr(2)
        moment2sum_thres =fltarr(2)
        moment3sum_thres =fltarr(2)
        expectedx = fltarr(num_links+1)
        expectedy = fltarr(num_links+1)
        x_link_offset=intarr(num_links+1)
        y_link_offset=intarr(num_links+1)
        linked_xerror = fltarr(num_links+1)
        linked_yerror = fltarr(num_links+1)
        
;        for k=0,num_subpixel-1 do begin
;            for l=0,num_subpixel-1 do begin
        for repeat_count=1,num_repeats do begin
            if (keyword_set(random_psf)) then begin
                ;We need to create a PSF at a random location in the shutter
                xrandom=10*randomu(seed,/uniform)
                yrandom=10*randomu(seed,/uniform)
;                print,xrandom,yrandom
                filename='interpolated/'+strtrim(string(repeat_count),2)+'.fits'
                print,filename
                writefits,filename,inimage
            endif
            for link_count=0,num_links do begin ; number of linked targets

                k=round(20*randomu(seed,/uniform)-0.5)
;       xoff=22 when the star is in the center of the pixel
;                k=10
                xoff = k+12 
                expectedx[link_count] = binned_peak[0] - 1.0*(xoff-2)/scale_factor - (num_subpixel/2)/scale_factor
;                print,expectedx[link_count],k,xoff
                l=round(20*randomu(seed,/uniform)-0.5)
;        yoff=20 when the star is in the center of the pixel
;                l=0
                yoff = l+10
;                print,' k ',k,' l ',l,' xoff ',xoff,' yoff',yoff
                expectedy[link_count] = binned_peak[1] - 1.0*yoff/scale_factor -(num_subpixel/2)/scale_factor
;                print,expectedy[link_count],l,yoff
                makepixelatedimage, inimage, outimage, total_counts, scale_factor, xoff, yoff, read_noise, seed, poisson_noise=poisson_noise
 ;               filename='miri_binned_'+strtrim(string(k),2)+'_'+strtrim(string(l),2)+'.fits'
;				print,filename
;				writefits,filename,outimage
                outimage = outimage + background
                find_peak, outimage, centroidbox,window_size,x_peak,y_peak
                find_peak, outimage, varcentroidbox,window_size+1,x_peak,y_peak
                add_bad_pixels,centroidbox,bad_pixel_probability,bad_value,seed,pixel_mask
                                ;	print,centroidbox
                if window_size eq 2 then begin
                    gwxbias =  -0.00826744
                    gwybias = +0.0154722
                    regxbias =-0.0125828
                    regybias =+0.0189102
                    corxbias =-0.0050545
                    corybias =+0.0189102
                    varxbias =-0.00563732
                    Varybias = 0.0128993
                    correctionfactor=0.472
                endif else begin
                    gwxbias =  -0.00610715
                    gwybias =  +0.0134574
                    corxbias = -0.00832816
                    corybias = +0.0153052
                    regxbias = -0.00832816
                    regybias = +0.0153052
                    varxbias = -0.00596658
                    varybias =  +0.0121879
                    correctionfactor=0.266
                endelse

                find_gw_moment, centroidbox, threshold, moment1,moment2,moment3,skewnesss,fwhm,iter_thres
                
                xcentroid = x_peak - (window_size) +moment1[0] +gwxbias
                ycentroid = y_peak - (window_size) +moment1[1] +gwybias
                gwxerror = expectedx-xcentroid
                gwyerror = expectedy-ycentroid
                find_moment, centroidbox, threshold, moment1,moment2,moment3,skewnesss
                xpixel_offset = moment1[0] - round(moment1[0])
                ypixel_offset = moment1[1] - round(moment1[1])
                xcorrection = correctionfactor * xpixel_offset
                ycorrection = correctionfactor * ypixel_offset
                xcentroid_cor = x_peak - (window_size) +moment1[0] +corxbias + xcorrection
                ycentroid_cor = y_peak - (window_size) +moment1[1] +corybias + ycorrection
                xerror_cor = expectedx-xcentroid_cor
                yerror_cor = expectedy-ycentroid_cor
                xcentroid = x_peak - (window_size) +moment1[0] +regxbias
                ycentroid = y_peak - (window_size) +moment1[1] +regybias
                xerror = expectedx-xcentroid
                yerror = expectedy-ycentroid
                find_varwindow_moment, varcentroidbox, threshold, moment1,moment2,moment3,skewnesss,iter_thres
                xcentroid = x_peak - (window_size+1) +moment1[0]  +varxbias
                ycentroid = y_peak - (window_size+1) +moment1[1]  +varybias
                varxerror = expectedx-xcentroid
                varyerror = expectedy-ycentroid
;                print,'raw'
            endfor              ;end of linked targets loop



;                   if (abs(xerror) gt 1.0) then begin
;                       print,expectedx,xcentroid,x_peak,expectedy,ycentroid,y_peak
;                        print,expectedx,xcentroid,expectedy,ycentroid,xerror,yerror,moment1,x_peak,y_peak,xoff,yoff

;                        print,expectedx,fitxloc,xcentroid,expectedy,fityloc,ycentroid,fitx,fity,xoff,yoff,fitxerror,fityerror
;                    endif
;            print,'moment2 ',moment2,' moment3',moment3
            ymodelfit = 0
            xmodelfit = 0
;            print,xmodelfit,ymodelfit
;            print,xcentroid,ycentroid
            printf,lun_out,expectedx[0],expectedy[0],xerror,yerror,gwxerror,gwyerror,varxerror,varyerror,xerror_cor,yerror_cor,FORMAT='(10F11.5)'
            xsum = xerror + xsum
            ysum = yerror + xsum
            varxsum = xerror + xsum
            varysum = yerror + xsum
            gwxsum = xerror + xsum
            gwysum = yerror + xsum
            x2sum = x2sum + xerror*xerror
            y2sum = y2sum + yerror*yerror
            x2sum_cor = x2sum_cor + xerror_cor*xerror_cor
            y2sum_cor = y2sum_cor + yerror_cor*yerror_cor
            gwx2sum = gwx2sum + gwxerror*gwxerror
            gwy2sum = gwy2sum + gwyerror*gwyerror
            varx2sum = varx2sum + varxerror*varxerror
            vary2sum = vary2sum + varyerror*varyerror
            count = count +1
        endfor
        rmsx = sqrt(x2sum/(count-1)) ; rms is about zero mean since 0 is the correct answer
        rmsy = sqrt(y2sum/(count-1))
        rmsxcor = sqrt(x2sum_cor/(count-1)) ; rms is about zero mean since 0 is the correct answer
        rmsycor = sqrt(y2sum_cor/(count-1))
        gwrmsx = sqrt(gwx2sum/(count-1)) ; rms is about zero mean since 0 is the correct answer
        gwrmsy = sqrt(gwy2sum/(count-1))
        varrmsx = sqrt(varx2sum/(count-1)) ; rms is about zero mean since 0 is the correct answer
        varrmsy = sqrt(vary2sum/(count-1))
        avgx = xsum/count
        avgy = ysum/count
        avgmoment2 = moment2sum/count
        avgmoment3 = moment3sum/count
        print,'1st moment x rms',rmsx
        print,'1st moment y rms',rmsy
        print,'gw moment x rms',gwrmsx
        print,'gw moment y rms',gwrmsy
        print,'var win moment x rms',varrmsx
        print,'var win moment y rms',varrmsy
        print,'corrected 1st mom x rms',rmsxcor
        print,'corrected 1st mom y rms',rmsycor
print,'total counts ',total_counts
;print,'avg xrms moment with threshold ',total(xrms_thres[where(xrms_thres ne 0)])/(2*(endx-startx+1)*(endy-starty+1)),' avg yrms moment with threshold',total(yrms_thres[where(yrms_thres ne 0)])/(2*(endx-startx+1)*(endy-starty+1))
;print,'avg xrms moment ',total(xrms)/(2*(endx-startx+1)*(endy-starty+1)),' avg yrms moment ',total(yrms)/(2*(endx-startx+1)*(endy-starty+1))


close,lun_out
free_lun,lun_out
end
