
;SHOW ME THE MATH:
;See the 
;https://courses.engr.illinois.edu/ece417/fa2020/slides/lec03.pdf
;
;Code reference:
;https://www.weatherclasses.com/uploads/1/3/1/3/131359169/fourier_analysis_using_idl.pdf
;

pro ex1_whitenoise_fft

repeat_n=10000
;repeat the experiment for [repeat_n] times
data_n=500
;number of data points
sample_arr=fltarr(data_n,repeat_n)
scale_n=500
amp_arr=fltarr(scale_n,repeat_n) & pha_arr=amp_arr
power_arr=fltarr(scale_n/2+1,repeat_n)
dt=1.
lag1=0.
sample_arr=randomn(9,data_n,repeat_n)*3
for i=0,repeat_n-1 do begin
  seedi=i
  ;sample_arr[*,i]=randomn(seedi,data_n)
  ft=fft(sample_arr[*,i])
  amp_arr[*,i]=abs(ft)
  
  power_arr[*,i]=(abs(ft[0:scale_n/2]))^2
  power_arr[0:scale_n/2-1,i]=power_arr[0:scale_n/2-1,i]*2.0
  ;why "scale_n/2"? See 
  ;https://www.weatherclasses.com/uploads/1/3/1/3/131359169/fourier_analysis_using_idl.pdf
  ;The remaining elements, up to index N/2 contain the power in
  ;the positive frequencies in ascending order. Past N/2, the elements contain the negative
  ;frequencies in descending order. In practice, I convert the spectrum to power using only the
  ;positive frequencies:
  ;IDL> power = (ABS(spectrum[0:N/2]))^2
  ;And, since I typically deal with FFTs when I’m doing regression analysis, I convert this
  ;power to variance explained at each frequency. Note that I do not multiply the Nyquist by
  ;2:
  ;IDL> power[0:N/2-1] = 2.0*power[0:N/2-1]
  pha_arr[*,i]=atan(ft,/phase)
endfor

; data_n is an integer giving the number of elements in a particular dimension
; dT is a floating-point number giving the sampling interval
ind = FINDGEN((data_n - 1)/2) + 1
dT = 1./data_n
is_N_even = (data_n MOD 2) EQ 0
if (is_N_even) then $
  freq = [0.0, ind, data_n/2, -data_n/2 + ind]/(data_n*dT) $
else $
  freq = [0.0, ind, -(data_n/2 + 1) + ind]/(data_n*dT)

;Show the amplitude spectrum [amp_spectral] vs [freq]
amp_spectral=mean(amp_arr,dim=2)
plot,freq,amp_spectral,psym=3
;Note that the figure is symetric with frequency=0
print,'Total of mean amplitude spectrum (meaningless): ',mean(amp_spectral)^2*scale_n

;Show the amplitude spectrum [amp_spectral^2] vs [freq]
ampsqr_spectral=mean(amp_arr^2,dim=2)
plot,freq,ampsqr_spectral,psym=3
;Note that the figure is symetric with frequency=0
print,'Total of mean amplitude^2 spectrum (supposed to be 1): ',mean(ampsqr_spectral)*scale_n

;Show the amplitude spectrum [amp_spectral] vs [freq] with only 
power_spectral=mean(power_arr,dim=2)
plot,freq(0:scale_n/2),power_spectral,psym=3, charsize=1.5,$
  title='The Amp_sqr (power) of a white noise case',$
  xtitle='Freq. index (k)',ytitle='Amp_sqr (i.e., power)'
print,'Total of mean power spectrum: ',mean(power_spectral)*(scale_n/2)

;After normalizing, the discrete Fourier power spextrum is Pk = TC_eq16
;Without normalizing, Pk = Ak^2 = [ TC_eq16 ] * (2 * sigma^2 / N)
;For white noise, TC_eq16=1
predict_aksqr = 2.*stddev(sample_arr)^2/data_n
oplot,!X.crange,predict_aksqr+[0,0]
pause


END

