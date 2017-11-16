pro find_peak,bigimage,centroidbox,window_size,x_peak,y_peak
insize = size(bigimage)
bigsize=insize[1]
hsize = fix(window_size/2.)
tbigimage = median(bigimage,window_size,/even)
;fits_write,'miri_med_image.fits',tbigimage
max_value=max(tbigimage[hsize:bigsize-hsize-1,hsize:bigsize-hsize-1],I)
ix = I MOD (bigsize - 2*hsize)
iy = I/(bigsize - 2*hsize)
centroidbox=fltarr(window_size,window_size)
if (ix gt window_size and ix+window_size lt bigsize and iy gt window_size and iy+window_size lt bigsize) then begin
    centroidbox=bigimage[ix-window_size:ix+window_size,iy-window_size:iy+window_size]
    x_peak=ix + hsize
    y_peak=iy + hsize
endif else begin
    centroidbox[*,*]=1
    x_peak=0
    y_peak=0
endelse

end
