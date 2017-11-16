
function get_miri_pixelated_image,fwhm,shift_x,shift_y,subsample,n_full,exp_center, $
  write_images=write_images

; create PSF from gaussian
if (not keyword_set(use_fits)) then begin
    n_sub_full = n_full*subsample
    center = fix(n_sub_full/2.*[1,1] + [shift_x,shift_y]*subsample)
;    print,center
    subsamp_image = psf_gaussian(npixel=n_sub_full,fwhm=fwhm*subsample,centroid=center)

    if (keyword_set(write_images)) then fits_write,'shifted_orig.fits',subsamp_image

    image = rebin(subsamp_image,n_full,n_full)

    exp_center = center/float(subsample)
endif else begin ; read in PSF image
    fits_read,'mips_24_3000K.astrom.fits',image,header
    image_size = size(image)
    getrot,header,rot,cdelt
    pixscale = abs(cdelt[0])*3600.

; shift the image if asked
    if ((n_elements(shift_x) GT 0L) AND (n_elements(shift_y) GT 0L)) then begin
        use_shift_x = fix(shift_x/((fwhm/6.)/6.))
        use_shift_y = fix(shift_y/((fwhm/6.)/6.))
;    print,'use = ', use_shift_x, use_shift_y

        nimage = replicate(min(image),2*image_size[1],2*image_size[2])
        tx = image_size[1]/2
        ty = image_size[2]/2
        nimage[tx+use_shift_x:tx+image_size[1]-1,ty+use_shift_y:ty+image_size[2]-1] = $
          image[0:image_size[1]-1-use_shift_x,0:image_size[2]-1-use_shift_y]
        image = nimage
        image_size = size(image)
        sxaddpar,header,'NAXIS1',image_size[1]
        sxaddpar,header,'NAXIS2',image_size[2]
    
        if (keyword_set(write_images)) then fits_write,'shifted_orig.fits',image
    endif

; rebin image so the PSF is the right size/pixel
;  the extra factor of 6 is to get 0.1 pixels
;  original pixelscale is ~0.5"/pixel
    hrebin,image,header,outsize=[image_size[1]*(fwhm/6.)/6.,image_size[2]*(fwhm/6.)/6.]
;print,'orig pixscale = ',pixscale
;print,'new pixscale = ',pixscale*(fwhm/6.)/6.  ; not correct
    image_size = size(image)

; make image bigger and fill with zeros
    if ((image_size[1] mod 2) EQ 0) then ext_off = 1 else ext_off = 0
    full_image = replicate(min(image),n_full,n_full)
    x1 = n_full/2 - image_size[1]/2 + ext_off + off_x
    x2 = n_full/2 + image_size[1]/2 + off_x
    y1 = n_full/2 - image_size[2]/2 + ext_off + off_y
    y2 = n_full/2 + image_size[2]/2 + off_y
    full_image[x1:x2,y1:y2] = image
    
    image = full_image
    image_size = size(image)
endelse

return, image

end
