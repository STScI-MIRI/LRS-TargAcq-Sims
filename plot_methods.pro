; program to illustrate the different methods of determining slopes on
; JWST

pro plot_methods,eps=eps

n_rej = 1
;total_eff = 0.25
exp_time = 930.
reads_per_sec = 1./30.
gain = 5.96
read_noise = 19. ; per 3s read
;read_noise = read_noise/sqrt(8) ; per 30s read (average of 8 reads)
;read_noise /= sqrt(12./34.) ; correcting to a single 30 sec read
dark_current = 0.3 ; e/sec
;sky_bkg = 50. ; photons/sec/micron/pixel
;sky_bkg *= total_eff  ; accounting for system throughput/QE

cr_prob_per_sec_per_pixel = 0.5/1000.       ;  50% of pixels affected in 1000s
cr_prob_per_int_per_pixel = cr_prob_per_sec_per_pixel*exp_time

sky_bkg = 6.76                  ; in e/sec
sky_bkg_out = (sky_bkg/reads_per_sec)/gain

create_ramp,200.,10.,exp_time,reads,dn_vals,dn_unc_vals,read_noise=read_noise,gain=gain, $
            reads_per_sec=reads_per_sec,seed=seed
n_dn_vals = n_elements(dn_vals)

dn_vals /= 1e3

; add a cosmic ray
cr_amp = 40000./gain            ; in DN
cr_amp /= 1e3
read_num = 18
dn_vals[read_num:n_dn_vals-1] += cr_amp

ps_file = 'miri_slope_methods'
setup_ps_output,ps_file,eps=eps

setup_colors,base_color,back_color,blue_color,red_color,green_color, $
             yel_color,purple_color,light_blue_color,grey_color, $
             line_color=line_color,red=red,green=green,blue=blue, $
             bw=bw

xrange = krange(reads/reads_per_sec)
yrange = krange(dn_vals)

kplot,[1],[1],/nodata,psym=1,xrange=xrange,yrange=yrange, $
      xtitle='time [sec]',ytitle='DN/10!U3!N',color=base_color,background=back_color
koplot,reads/reads_per_sec,dn_vals,psym=1,color=base_color

dn_vals_rebin = rebin(dn_vals,4)
reads_rebin = rebin(reads,4)
slope_vals = ([dn_vals_rebin[2:3]] - [dn_vals_rebin[1:2]])/((n_dn_vals/4.)*[1.,2.])

koplot,reads_rebin[1:3]/reads_per_sec,dn_vals_rebin[1:3],psym=3,symsize=3.,color=red_color
koplot,reads_rebin[1:3]/reads_per_sec,dn_vals_rebin[1:3],psym=100,linestyle=2,color=red_color,thick=4

bindxs = (indgen(4))*fix((n_dn_vals/4.))
eindxs = (indgen(4) + 1)*fix((n_dn_vals/4.)) - 1
slope_vals = (dn_vals[eindxs] - dn_vals[bindxs])/(n_dn_vals/4.)

tindxs = [0,read_num-1]
koplot,reads[tindxs]/reads_per_sec,dn_vals[tindxs],psym=100,linestyle=0,color=green_color,thick=4

tindxs = [bindxs[1],eindxs[1]]
koplot,reads[tindxs]/reads_per_sec,dn_vals[tindxs],psym=100,linestyle=1,color=blue_color,thick=4
tindxs = [bindxs[2],eindxs[2]]
koplot,reads[tindxs]/reads_per_sec,dn_vals[tindxs],psym=100,linestyle=1,color=blue_color,thick=4
tindxs = [bindxs[3],eindxs[3]]
koplot,reads[tindxs]/reads_per_sec,dn_vals[tindxs],psym=100,linestyle=1,color=blue_color,thick=4

close_ps_output,eps=eps

end
