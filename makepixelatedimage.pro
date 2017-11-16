Pro makepixelatedimage, inimage, outimage, total_counts, scale_factor, x_offset, y_offset, read_noise, seed, poisson_noise=poisson_noise

insize = size(inimage)
inx = insize[1]
iny = insize[2]

outx = floor(inx / scale_factor +0.5)
outy = floor(iny / scale_factor +0.5)

outimage=fltarr(outx, outy)

  if (total(~finite(inimage)) ne 0) then begin
      inimage[where(~finite(inimage))]=0.0
      inimage[where(inimage lt -1e6)]=0
      print,total(inimage)
  endif
scaledin = inimage/total(inimage)

for i=0,outx-3 do begin
    for j=0,outy-3 do begin
        expected_counts = total(scaledin[x_offset+i*scale_factor:x_offset+(i+1)*scale_factor,y_offset+j*scale_factor:y_offset+(j+1)*scale_factor])
        expected_counts = total_counts*expected_counts
        if (not keyword_set(poisson_noise)) then begin
            actual_counts = expected_counts
        end else begin
;            actual_counts = poidev(expected_counts,seed=seed)
            actual_counts = randomu(seed,poisson=expected_counts)
;            if (actual_counts gt 1) then begin
;                print,'adding poisson noise',expected_counts,actual_counts,i,j
;            endif
        endelse
                                ;add readnoise
        actual_counts = actual_counts + read_noise*Randomu(seed,/normal)
;        if (expected_counts gt 1) then begin
;            print,'adding read noise',expected_counts,actual_counts
;        endif
        outimage[i,j]=actual_counts
    endfor
endfor
if (total((~FINITE(outimage))) ne 0) then begin
   www=xxx
endif
;print,'output of make pixletated',total(outimage)
end

