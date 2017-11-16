
pro plot_miri_runs,eps=eps

files = ['miri_run_h4_255_m2_r3.dat', $
         'miri_run_h3_15_m2_r2.dat', $
         'miri_run_h2_10_m2_r1.dat', $
         'miri_run_h2_5.6_m2_r1.dat']
n_files = n_elements(files)

for i = 0,(n_files-1) do begin
    plot_many_centroids,files[i],eps=eps,/radial
    if (not keyword_set(eps)) then begin
        ans = ''
        read,'Continue: ',ans
    endif

    plot_many_centroids,repstr(files[i],'.dat','_hs.dat'),eps=eps,/radial
    if (not keyword_set(eps)) then begin
        ans = ''
        read,'Continue: ',ans
    endif
endfor

end
