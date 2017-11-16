Pro find_varwindow_moment,centroidbox,iter_thres,max_iter,moment1,moment2,moment3,skewness
Insize = size(centroidbox)
winsize=insize[1]
boxsize=(insize[1]-3)/2
moment1=fltarr(2)
moment2=fltarr(2)
moment3=fltarr(2)
skewness=fltarr(2)
rawmoment2=fltarr(2)
rawmoment3=fltarr(2)

sum = 0.
xsum=0.
x2sum=0.
x3sum=0.
ysum=0.
y2sum=0.
y3sum=0.
;Make a first guess at where the peak is using an unweighted moment
for i=1,winsize do begin
    for j=1,winsize do begin
        xsum=i*centroidbox[i-1,j-1] +xsum
        x2sum=i*i*centroidbox[i-1,j-1] +x2sum
        x3sum=i*i*i*centroidbox[i-1,j-1] +x3sum
        ysum=j*centroidbox[i-1,j-1] +ysum
        y2sum=j*j*centroidbox[i-1,j-1] +y2sum
        y3sum=j*j*j*centroidbox[i-1,j-1] +y3sum
    endfor
endfor
moment1[0]=xsum/total(centroidbox) - 1
moment1[1]=ysum/total(centroidbox) - 1
;print,moment1[0],moment1[1]
oldmoment1=[0,0]
;Now iterate the gaussian moment until the peak doesn't change by some threshold
count=0
While ((abs(moment1[0]-oldmoment1[0]) gt iter_thres) or (abs(moment1[1]-oldmoment1[1] gt iter_thres))and count lt max_iter) do begin
;    print,count,moment1
    sum=0
    xsum=0.
    x2sum=0.
    x3sum=0.
    ysum=0.
    y2sum=0.
    y3sum=0.
    sumweight=0
    for i=1,winsize do begin
        for j=1,winsize do begin
            xweight=0
            yweight=0
            xoff = (i-1)-moment1[0]
            yoff = (j-1)-moment1[1]
            if (abs(xoff)) le (boxsize) then begin
                xweight=1
            endif else begin
                                ;this will be for pixels that need to be downweighted, they would have been in the original window
                if (abs(xoff) gt boxsize) and (abs(xoff) lt  boxsize+1 ) then begin
                    xweight = boxsize+1 - abs(xoff)
                endif 
            endelse
            
            if (abs(yoff)) le (boxsize) then begin
                yweight=1
            endif else begin
                                ;this will be for pixels that need to be downweighted, they would have been in the original window
                if (abs(yoff) gt boxsize) and (abs(yoff) lt  boxsize+1 ) then begin
                    yweight = boxsize+1 - abs(yoff)
                endif 
            endelse
            
            weight=xweight*yweight
;            print,i,j,weight,xweight,yweight,xoff,yoff,moment1[0],moment1[1]
            sumweight=sumweight + weight
            sum=centroidbox[i-1,j-1]*weight + sum
            xsum=i*centroidbox[i-1,j-1]*weight +xsum
            x2sum=i*i*centroidbox[i-1,j-1]*weight +x2sum
            x3sum=i*i*i*centroidbox[i-1,j-1]*weight +x3sum
            ysum=j*centroidbox[i-1,j-1]*weight +ysum
            y2sum=j*j*centroidbox[i-1,j-1]*weight +y2sum
            y3sum=j*j*j*centroidbox[i-1,j-1]*weight +y3sum
;            print,i,j,sum,xsum,ysum,weight
        endfor
    endfor
    oldmoment1=moment1
    moment1[0]=xsum/sum - 1
    moment1[1]=ysum/sum - 1
    count=count+1
endwhile
;print,moment1[0],moment1[1],sumweight,count,boxsize

;calculate the 1st, 2nd, and 3rd raw moments
rawmoment2[0]=x2sum/sum
rawmoment2[1]=y2sum/sum
rawmoment3[0]=x3sum/sum
rawmoment3[1]=y3sum/sum
; We are interested in the 2nd and 3rd central moments
; From Wolfram's web page they are defined in terms of the raw moments
moment2[0]=rawmoment2[0]-moment1[0]*moment1[0]
moment2[1]=rawmoment2[1]-moment1[1]*moment1[1]
moment3[0]= 2 * moment1[0]^3 - 3*moment1[0]*rawmoment2[0] + rawmoment3[0]
moment3[1]= 2 * moment1[1]^3 - 3*moment1[1]*rawmoment2[1] + rawmoment3[1]
;The skewness depends on both the 2nd and 3rd central moments
skewness=moment3/moment2^1.5
end
