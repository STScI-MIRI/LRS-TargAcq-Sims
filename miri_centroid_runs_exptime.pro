; all the runs needed for the technical report

pro miri_centroid_runs_exptime,n_test,wavelength,filter,req_val=req_val

if (not keyword_set(req_val)) then req_val = 1.
if (filter EQ 'ND') then nd = 1 else nd = 0
test_many_centroids,wavelength,n_test,method=2,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime960.dat',exp_time=960.,nd=nd
test_many_centroids,wavelength,n_test,method=2,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime720.dat',exp_time=720.,nd=nd
test_many_centroids,wavelength,n_test,method=2,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime480.dat',exp_time=480.,nd=nd
test_many_centroids,wavelength,n_test,method=2,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime360.dat',exp_time=360.,nd=nd
test_many_centroids,wavelength,n_test,method=2,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime240.dat',exp_time=240.,nd=nd
test_many_centroids,wavelength,n_test,method=3,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime120.dat',exp_time=120.,nd=nd
test_many_centroids,wavelength,n_test,method=3,req_val=1.,filename='miri_run_h2_'+filter+'_m2_r1_exptime12.dat',exp_time=12., $
                    reads_per_sec=1./3.,nd=nd

end
