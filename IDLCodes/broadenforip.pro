; NAME;
;   broadenforip.pro
;
; AUTHOR:
;   Hiroyuki Tako Ishikawa (contact: hishikaw@uwo.ca)
;
; PURPOSE:
;   To broaden spectra with instrumental profiles (IP)
;
; INPUT PARAMETERS:
;   wavelength (µm)
;   flux
;   
; RETURNS:
;   wavelengths convolved with IP
;
; CALLING SEQUENCE:
;   flux_broadened = broadenforip(wavelength, flux)
;

function inverted_gaussian_with_trend, x, amp, cen, wid, offset, slope
  compile_opt idl2
  return, amp * exp(-(x - cen) ^ 2 / (2 * wid ^ 2)) + offset + slope * x
end


function convolve_spectrum, w_all, f_all, lsf_sigma_in_micron, length_lsf_to_convolve
  compile_opt idl2
  if n_elements(length_lsf_to_convolve) eq 0 then length_lsf_to_convolve = 100

  ; lsf_sigmaの計算
  lambda_per_pixel_model = (max(w_all) - min(w_all)) / (n_elements(w_all) - 1)
  N_pixel_model = length_lsf_to_convolve - 1
  lsf_sigma = lsf_sigma_in_micron / (lambda_per_pixel_model * N_pixel_model)

  x_lsf_full = findgen(length_lsf_to_convolve) / (length_lsf_to_convolve - 1)
  lsf_full = inverted_gaussian_with_trend(x_lsf_full, 1, 0.5, lsf_sigma, 0, 0)
  lsf_full = lsf_full / total(lsf_full)
  f_all_convolved = convol(double(f_all), double(lsf_full), /edge_truncate)
  return, f_all_convolved
end


function broadenforip, wavelength, flux, wav_micron=wav_micron
  compile_opt idl2

  if n_elements(wav_micron) eq 0 then wav_micron = mean(wavelength)  
  sigma_ip = get_sigma_ip(wav_micron)

  return, convolve_spectrum(wavelength, flux, sigma_ip)
end