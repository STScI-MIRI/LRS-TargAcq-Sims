
pro plot_many_centroids,filename,eps=eps,radial=radial,req_val=req_val, $
  est_req_flux=est_req_flux

setup_colors,base_color,back_color,blue_color,red_color,green_color, $
             yel_color,purple_color,light_blue_color, $
             line_color=line_color,red=red,green=green,blue=blue, $
             bw=bw

pixscale = 0.11

if (n_elements(filename) EQ 0L) then filename = 'miri_test_centroids.dat'

ps_file = repstr(filename,'.dat','_plot')
if (keyword_set(radial)) then ps_file += '_radial'
setup_ps_output,ps_file,eps=eps

readcol,filename,fluxes,c1,c1u,c2,c2u,ng,ec1,ec1u,ec2,ec2u,eng
c1u *= 1e3
c2u *= 1e3
ec1u *= 1e3
ec2u *= 1e3
fluxes *= 1e3

xrange = krange(fluxes,kplot_type='o')
yrange = krange([c1u,c2u,ec1u,ec2u]*pixscale,kplot_type='o')
yrange = [0.0001,0.1]*1e3

kplot,[1],[1],/no_data,xrange=xrange,yrange=yrange,ystyle=8, $
      xtitle='flux [!4l!3Jy]',ytitle='uncertainty [mas]', $
      color=base_color,background=back_color,kplot_type='oo', $
      position=[0.15,0.15,0.9,0.95]

if (not keyword_set(req_val)) then req_val = 0.005*1e3
if (keyword_set(radial)) then begin
    rad_slope = sqrt(c1u^2 + c2u^2)
    rad_est_slope = sqrt(ec1u^2 + ec2u^2)

    ; determine where the results cross the requirements
    indxs = where((rad_est_slope*pixscale LT req_val*5.) AND (rad_est_slope*pixscale GT req_val*0.2),n_indxs)
    est_req_flux = 0.0
    if (n_indxs GT 0) then begin
        print,'required centroiding accuracy = ',req_val
        req_val_full = interpol(fluxes[indxs],rad_slope[indxs]*pixscale,[req_val])
        print,'full_req_flux = ',req_val_full
        req_val_est = interpol(fluxes[indxs],rad_est_slope[indxs]*pixscale,[req_val])
        print,'est_req_flux = ',req_val_est
        est_req_flux = req_val_est
    endif

    koplot,fluxes,rad_slope*pixscale,psym=2,color=blue_color
    koplot,fluxes,rad_est_slope*pixscale,psym=4,color=red_color
    
    ktype = [[2,2],[0,0],[4,4],[0,0]]
    klabel = ['fit','fit fraction','est', $
              'est fraction']
    line_color = [blue_color,blue_color,red_color,red_color]

;    ktype2 = [[-1,0],[-2,0]]
;    klabel2 = ['requirement','centroiding']
;    line_color2 = [green_color,green_color]
;    klegend,[0.6,0.7],ktype2,klabel2,color=base_color,kplot_type='oo',charsize=1.3,/box, $
;            line_color=line_color2
endif else begin
    koplot,fluxes,c1u*pixscale,psym=1,color=blue_color
    koplot,fluxes,c2u*pixscale,psym=2,color=blue_color
    koplot,fluxes,ec1u*pixscale,psym=3,color=red_color
    koplot,fluxes,ec2u*pixscale,psym=4,color=red_color
    
    ktype = [[1,1],[2,2],[3,3],[4,4]]
    klabel = ['x (slope fit)','y (slope fit)','x (slope est)','y (slope est)']
    line_color = [blue_color,blue_color,red_color,red_color]
endelse

;koplot,xrange,*1e3*[1.,1.],color=green_color,linestyle=1
koplot,xrange,req_val*[1.,1.],color=green_color,linestyle=2

klegend,[0.2,0.3],ktype,klabel,color=base_color,kplot_type='oo',charsize=1.3,/box, $
        line_color=line_color

; now plot the 2nd y-axis and data

;AXIS, YAXIS=1, YLOG=1, YRANGE=[0.1, 100], /SAVE
axis,yaxis=1,ylog=0,yrange=[0.,1.05],/save,color=base_color,ystyle=1, $
     ythick=2.,charthick=2.,ytitle='Fraction Successful',charsize=1.5
koplot,fluxes,ng/max(ng),color=blue_color,linestyle=0,psym=100
koplot,fluxes,eng/max(eng),color=red_color,linestyle=0,psym=100
;koplot,fluxes,ng*0.09*1e3/max(ng),color=blue_color,linestyle=0,psym=100
;koplot,fluxes,eng*0.09*1e3/max(eng),color=red_color,linestyle=0,psym=100

close_ps_output,eps=eps

end
