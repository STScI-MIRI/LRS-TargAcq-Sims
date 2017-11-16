; all the runs needed for the technical report

pro miri_centroid_runs,n_test

;test_many_centroids,5.6,n_test,method=2,req_val=1.,filename='miri_run_h2_5.6_m2_r1.dat'
;test_many_centroids,5.6,n_test,method=2,req_val=1.,filename='miri_run_h2_5.6_m2_r1_hs.dat',/highsky

;test_many_centroids,25.5,n_test,method=2,req_val=3.,filename='miri_run_h4_255_m2_r3_hs.dat',highsky=1
;test_many_centroids,25.5,n_test,method=2,req_val=3.,filename='miri_run_h4_255_m2_r3.dat'

test_many_centroids,15.0,n_test,method=2,req_val=2.,filename='miri_run_h3_15_m2_r2_hs.dat',/highsky
test_many_centroids,10.0,n_test,method=2,req_val=1.,filename='miri_run_h2_10_m2_r1_hs.dat',/highsky

test_many_centroids,15.0,n_test,method=2,req_val=2.,filename='miri_run_h3_15_m2_r2.dat'
test_many_centroids,10.0,n_test,method=2,req_val=1.,filename='miri_run_h2_10_m2_r1.dat'

end
