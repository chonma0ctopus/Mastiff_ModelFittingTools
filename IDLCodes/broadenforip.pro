; NAME;
; broadenforip.pro
;
; AUTHOR:
; Hiroyuki Tako Ishikawa (contact: hishikaw@uwo.ca)
;
; PURPOSE:
; To broaden spectra with IGRINS's instrumental profile
;
; INPUT PARAMETERS:
;
; RETURNS:
;
; CALLING SEQUENCE:
;
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
  lsf_full = lsf_full / total(lsf_full) ; 正規化
  f_all_convolved = convol(double(f_all), double(lsf_full), /edge_truncate)
  return, f_all_convolved
end


function broadenforip, wavelength, flux, wav_micron=wav_micron
  compile_opt idl2

  if n_elements(wav_micron) eq 0 then wav_micron = mean(wavelength)  
  sigma_ip = get_sigma_ip(wav_micron)

  return, convolve_spectrum(wavelength, flux, sigma_ip)
end

; test
; function testrun
;   compile_opt idl2

;   ; データの読み込み
;   file_path = '/Users/chonmac/Downloads/test_spec_broadenforip.dat'
;   readcol, file_path, wavelength, flux, format = 'D,D'

;   ; broadenforipを使ってブロードニング
;   flux_broadened = broadenforip(wavelength, flux)

;   ; グラフィックスデバイスを24ビットカラーモードに設定
;   device, decomposed = 1, true_color = 24

;   ; 結果のプロット（背景を白に設定）
;   PLOT, wavelength, flux, $
;     xtitle = 'Wavelength', $
;     ytitle = 'Flux', $
;     title = 'Original vs Broadened Spectrum', $
;     /nodata, $
;     background = 'ffffff'x, $ ; 背景色を白に設定 (16進数)
;     color = '000000'x ; テキストと軸の色を黒に設定 (16進数)

;   ; オリジナルのスペクトルをプロット
;   oplot, wavelength, flux, color = '000000'x ; 黒色 (16進数)

;   ; ブロードニングされたスペクトルをプロット
;   oplot, wavelength, flux_broadened, color = 'ff5f7f'x ; 赤色 (16進数)

;   ; legendの追加
;   AL_LEGEND, ['Original', 'Broadened'], $
;     linestyle = [0, 0], $
;     colors = ['000000'x, 'ff5f7f'x], $ ; 黒と赤 (16進数)
;     position = [0.1, 0.9], $
;     /normal, $
;     charsize = 1.0

;   ; plot, wavelength, flux, xtitle = 'Wavelength', ytitle = 'Flux', title = 'Original vs Broadened Spectrum'
;   ; oplot, wavelength, flux_broadened, color = 255 ; 赤色のインデックス
;   ; cgLegend, titles = ['Original', 'Broadened'], linestyle = [0, 0], colors = [0, 255], /addcmd
; end