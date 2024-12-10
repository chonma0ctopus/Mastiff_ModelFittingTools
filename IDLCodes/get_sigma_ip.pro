; NAME;
; get_sigma_ip.pro
;
; AUTHOR:
; Hiroyuki Tako Ishikawa (contact: hishikaw@uwo.ca)
;
; PURPOSE:
; To obtain sigma (µm) of instrumental broadening profiles (IP) for IGRINS
;
; INPUT PARAMETERS:
;    wav_micron: wavelengths (µm)
;    
; RETURNS:
;    Gaussian sigma of the instrumental profile (IP) broadening (µm)
;    Note that The current (Dec, 2024) default settings are optimized 
;        for the April-June 2018 Gemini South/IGRINS IP. Users will 
;        need to edit this file if they are using slow rotators' 
;        spectra from other instruments. (Ideally, it will be 
;        possible to adjust it from an external file, such as config 
;        files, in the future)
;
; CALLING SEQUENCE:
;    sigma_ip = get_sigma_ip(wav_micron)
;

function linear_func, x, a, b
  compile_opt idl2
  return, a * x + b
end

function get_sigma_ip, wav_micron
  compile_opt idl2
  ; wav_micron: wavelengths in the order centre (µm)

  switching_point = 1.850
  popt1 = [4.15512657e-05, -4.16661217e-05]
  popt2 = [4.25428496e-05, -5.86983801e-05]
  const_adjust_lowerlim_H = 0 ; 0.51e-5
  const_adjust_lowerlim_K = 0 ; 0.72e-5

  if wav_micron le switching_point then begin
    return, linear_func(wav_micron, popt1[0], popt1[1]) - const_adjust_lowerlim_H
  endif else begin
    return, linear_func(wav_micron, popt2[0], popt2[1]) - const_adjust_lowerlim_K
  endelse
end