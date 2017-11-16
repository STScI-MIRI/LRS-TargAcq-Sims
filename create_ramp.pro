
pro create_ramp,zero_point,slope,exp_time,reads,dn_vals,dn_unc_vals,dn_ran_unc_vals,dn_cor_unc_vals, $
                seed=seed,read_noise=read_noise,gain=gain, $
                random_only=random_only,cor_only=cor_only,reads_per_sec=reads_per_sec

; zero_point in DN, slope in 10^3 DN/sec

if (not keyword_set(gain)) then gain = 7.1
if (not keyword_set(read_noise)) then read_noise = 2.30*100.0  ; in electrons, sqrt(12/N) times worse than 8 sec value

if (not keyword_set(reads_per_sec)) then reads_per_sec = 8.0
n_reads = exp_time*reads_per_sec + 1
reads = findgen(n_reads)

rn_sigma = (read_noise/gain)*randomn(seed,n_reads,/normal,/double)

ave_rn = total(rn_sigma)/n_reads
std_rn = sqrt(total((rn_sigma - ave_rn)^2)/(n_reads-1))

photon_sigma = (sqrt(2.0*gain*slope/reads_per_sec)/gain)*randomn(seed,n_reads,/double,/normal)

dn_vals = dblarr(n_reads)
dn_vals[0] = zero_point
for i = 1,(n_reads-1) do begin
    dn_vals[i] = dn_vals[i-1] + slope/reads_per_sec
    if (not keyword_set(random_only)) then begin
        dn_vals[i] = dn_vals[i] + photon_sigma[i]
    endif
endfor

; add read noise
if (not keyword_set(cor_only)) then begin
    dn_vals = dn_vals + rn_sigma
endif

; get uncs
dn_unc_vals = dblarr(n_reads)
dn_ran_unc_vals = dblarr(n_reads)
dn_cor_unc_vals = dblarr(n_reads)
indxs = where(dn_vals GT 0.0,n_indxs)
if (n_indxs GT 0) then begin
    dn_ran_unc_vals[indxs] = read_noise/gain
    dn_cor_unc_vals[indxs] = sqrt(2.0*gain*dn_vals[indxs])/gain
    if (not keyword_set(random_only)) then begin
        dn_unc_vals[indxs] = sqrt((dn_cor_unc_vals)^2 + (dn_ran_unc_vals)^2)
    endif else begin
        dn_unc_vals[indxs] = dn_ran_unc_vals
    endelse
endif

indxs = where(dn_vals LE 0.0,n_indxs)
if (n_indxs GT 0) then begin
    dn_unc_vals[indxs] = read_noise/gain
endif

end
