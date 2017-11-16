
function kg_linfit_cor, x, y, in_unc, read_noise, yfit=yfit, sigma=sigma

gain = 7.1

n_x = n_elements(x)

;unc = sqrt(in_unc^2 + replicate(read_noise/gain,n_x)^2)
unc = in_unc

s = 0.0d0
sx = 0.0d0
sxx = 0.0d0
sy = 0.0d0
sxy = 0.0d0
for i = (n_x-1),1,-1 do begin
    sigma2 = 1.0/unc[i]^2
    s = s + sigma2
    sx = sx + x[i]*sigma2
    sxx = sxx + x[i]*x[i]*sigma2
    sy = sy + (y[i] - y[i-1])*s
    sxy = sxy + (y[i] - y[i-1])*sx
endfor
sigma2 = 1.0/unc[0]^2
s = s + sigma2
sx = sx + x[0]*sigma2
sxx = sxx + x[0]*x[0]*sigma2
sy = sy + y[0]*s
sxy = sxy + y[0]*sx

delt = s*sxx - sx^2
fit_coeff = fltarr(2)
fit_coeff[0] = (sxx*sy - sx*sxy)/delt
fit_coeff[1] = (s*sxy - sx*sy)/delt

yfit = fit_coeff[0] + fit_coeff[1]*x

; determine the coeff unc
sigma = fltarr(2)

new_s = 0.0d0
new_sx = 0.0d0
sigma_a = 0.0d0
sigma_b = 0.0d0
p_unc = sqrt(2.0*gain*fit_coeff[1])/gain
for i = (n_x-1),1,-1 do begin
    sigma2 = 1.0/unc[i]^2
    new_s = new_s + sigma2
    new_sx = new_sx + x[i]*sigma2

    sum_val = (sxx*new_s - sx*new_sx)^2
    sigma_a = sigma_a + p_unc^2*sum_val
    sum_val = (s*new_sx - sx*new_s)^2
    sigma_b = sigma_b + p_unc^2*sum_val
endfor
; y_1 term
;sigma_a = sigma_a + ((sxx - sx*x[0])/unc[0])^2
;sigma_b = sigma_b + ((s*x[0] - sx)/unc[0])^2
; final number
sigma_a = sqrt(sigma_a)/delt
sigma_b = sqrt(sigma_b)/delt

; compute read noise component
total_x = total(x)
total_x2 = total(x^2)
sigma_a_rn = total_x2*(read_noise/gain)^2/(n_x*total_x2 - total_x^2)
sigma_b_rn = n_x*(read_noise/gain)^2/(n_x*total_x2 - total_x^2)

sigma[0] = sqrt(sigma_a^2 + sigma_a_rn)
sigma[1] = sqrt(sigma_b^2 + sigma_b_rn)
;sigma[0] = sigma_a
;sigma[1] = sigma_b
;sigma[0] = sqrt(sigma_a_rn)
;sigma[1] = sqrt(sigma_b_rn)

return, fit_coeff

end
