;This code plot Fig 5 of Torrence and Compo, 1998 paper
;alfa=0 for fig 5a, alfa=0.7 for fig 5b
;Shuang Xu, 2022-05-08 

pro ex3_tc_fig5_wavelet

repeat_n=10000
alfa=.7 ;when (alfa EQ 0), white-noise; when (alfa NE 0), red-noise 
;repeat Monte Carlo test for repeat_n times (in the Torrence paper repeat_n=100000)
data_n=512 ;length of tiem series data
sample_arr=randomn(7,data_n,repeat_n) ;generate white noise
if (alfa NE 0) then begin
  for i=1,data_n-1 do begin
    sample_arr[i,*]=alfa*sample_arr[i-1,*]+sample_arr[i,*]
  endfor
endif

;define the scales according to smallest period [s0] and sub-octave resolution [dj]
dt=1.
s0 = dt*2.    ; this says start at a scale of dx
power1=0
dj = (1./2.)*(2^power1)  ; this will produce 8 sub-octaves per octave, 
;dj in Torrence and Compo, 1998 is equivalent to [1/J] in (5.13) in 
;Dr. Cao Chen's 2016 dissertation. It also denotes the resolution of spectral dimension
power2=fix(log2(data_n*dt/s0))
j1 = power2/dj  ; this says do 9 powers-of-two with dj sub-octaves each
scale_n=j1+1    
power_arr=fltarr(scale_n,repeat_n) & pha_arr=power_arr & signif_arr=power_arr & expect_arr=power_arr
power_arr95=fltarr(scale_n)

for i=0,repeat_n-1 do begin
  input_case=sample_arr[*,i]-mean(sample_arr[*,i]) ;make input with average 0.
  wave=wavelet_xu(input_case,dt,mother='morlet',scale=scale,S0=s0,lag1=alfa,$
    period=period,PAD=0,COI=coi,DJ=dj,J=j1,fft_theor=fft_theor,signif=signif,$
    xu_tceq7=1)
  ;signif = WAVE_SIGNIF(input_case,dt,scale,0,LAG1=alfa,MOTHER='morlet',fft_theor=fft_theor)
  signif_arr[*,i]=signif
  expect_arr[*,i]=fft_theor
    
  power=abs(wave)^2
  power_slice=power[fix(data_n/2),*]  ;data_n/2 = 256
  ;The local wavelet spectra were slices taken at time n = 256 out of N = 512 points.
  power_arr[*,i]=power_slice
  pha_arr[*,i]=atan(wave[fix(data_n/2),*],/phase)
endfor
for i=0,scale_n-1 do begin
  power_slice2=power_arr[i,*]
  power_slice2_sort=power_slice2(sort(power_slice2))
  power_arr95[i]=power_slice2_sort[repeat_n*.95]
endfor

stddev_sample=stddev(sample_arr);mean(stddev(sample_arr,dimension=1))
normalized_var=mean(power_arr,dim=2)/(stddev_sample^2);normalize the result
normalized_signif=mean(signif_arr,dimension=2)/(stddev_sample^2)
normalized_expect=mean(expect_arr,dimension=2)/(stddev_sample^2)

plot,period,normalized_var,/xlog,charsize=1.5,yrange=[0,4]*(alfa EQ 0?1:5),$
  xtitle='Period',ytitle='Normalized var',psym=4,$
  title='Fig 5a, Normalized mean-variance (repeat_n='+string(repeat_n,format='(E7.1)')+')'
oplot,period,power_arr95/(stddev(sample_arr)^2),psym=5

oplot,period,normalized_signif,linestyle=0
oplot,period,normalized_expect,linestyle=0
pause
END





;******************************************************************* WAVELET
;+
; NAME:   WAVELET
;
; PURPOSE:   Compute the WAVELET transform of a 1D time series.
;       
;
; CALLING SEQUENCE:
;
;      wave = WAVELET(Y,DT)
;
;
; INPUTS:
;
;    Y = the time series of length N.
;
;    DT = amount of time between each Y value, i.e. the sampling time.
;
;
; OUTPUTS:
;
;    WAVE is the WAVELET transform of Y. This is a complex array
;    of dimensions (N,J+1). FLOAT(WAVE) gives the WAVELET amplitude,
;    ATAN(IMAGINARY(WAVE),FLOAT(WAVE)) gives the WAVELET phase.
;    The WAVELET power spectrum is ABS(WAVE)^2.
;
;
; OPTIONAL KEYWORD INPUTS:
;
;    S0 = the smallest scale of the wavelet.  Default is 2*DT.
;
;    DJ = the spacing between discrete scales. Default is 0.125.
;         A smaller # will give better scale resolution, but be slower to plot.
;
;    J = the # of scales minus one. Scales range from S0 up to S0*2^(J*DJ),
;        to give a total of (J+1) scales. Default is J = (LOG2(N DT/S0))/DJ.
;
;    MOTHER = A string giving the mother wavelet to use.
;            Currently, 'Morlet','Paul','DOG' (derivative of Gaussian)
;            are available. Default is 'Morlet'.
;
;    PARAM = optional mother wavelet parameter.
;            For 'Morlet' this is k0 (wavenumber), default is 6.
;            For 'Paul' this is m (order), default is 4.
;            For 'DOG' this is m (m-th derivative), default is 2.
;
;    PAD = if set, then pad the time series with enough zeroes to get
;         N up to the next higher power of 2. This prevents wraparound
;         from the end of the time series to the beginning, and also
;         speeds up the FFT's used to do the wavelet transform.
;         This will not eliminate all edge effects (see COI below).
;
;    LAG1 = LAG 1 Autocorrelation, used for SIGNIF levels. Default is 0.0
;
;    SIGLVL = significance level to use. Default is 0.95
;
;    VERBOSE = if set, then print out info for each analyzed scale.
;
;    RECON = if set, then reconstruct the time series, and store in Y.
;            Note that this will destroy the original time series,
;            so be sure to input a dummy copy of Y.
;
;    FFT_THEOR = theoretical background spectrum as a function of
;                Fourier frequency. This will be smoothed by the
;                wavelet function and returned as a function of PERIOD.
;
;
; OPTIONAL KEYWORD OUTPUTS:
;
;    PERIOD = the vector of "Fourier" periods (in time units) that corresponds
;           to the SCALEs.
;
;    SCALE = the vector of scale indices, given by S0*2^(j*DJ), j=0...J
;            where J+1 is the total # of scales.
;
;    COI = if specified, then return the Cone-of-Influence, which is a vector
;        of N points that contains the maximum period of useful information
;        at that particular time.
;        Periods greater than this are subject to edge effects.
;        This can be used to plot COI lines on a contour plot by doing:
;            IDL>  CONTOUR,wavelet,time,period
;            IDL>  PLOTS,time,coi,NOCLIP=0
;
;    YPAD = returns the padded time series that was actually used in the
;         wavelet transform.
;
;    DAUGHTER = if initially set to 1, then return the daughter wavelets.
;         This is a complex array of the same size as WAVELET. At each scale
;         the daughter wavelet is located in the center of the array.
;
;    SIGNIF = output significance levels as a function of PERIOD
;
;    FFT_THEOR = output theoretical background spectrum (smoothed by the
;                wavelet function), as a function of PERIOD.
;
;
; [ Defunct INPUTS:
; [   OCT = the # of octaves to analyze over.           ]
; [         Largest scale will be S0*2^OCT.             ]
; [         Default is (LOG2(N) - 1).                   ]
; [   VOICE = # of voices in each octave. Default is 8. ]
; [          Higher # gives better scale resolution,    ]
; [          but is slower to plot.                     ]
; ]
;
;----------------------------------------------------------------------------
;
; EXAMPLE:
;
;    IDL> ntime = 256
;    IDL> y = RANDOMN(s,ntime)       ;*** create a random time series
;    IDL> dt = 0.25
;    IDL> time = FINDGEN(ntime)*dt   ;*** create the time index
;    IDL> 
;    IDL> wave = WAVELET(y,dt,PERIOD=period,COI=coi,/PAD,SIGNIF=signif)
;    IDL> nscale = N_ELEMENTS(period)
;    IDL> LOADCT,39
;    IDL> CONTOUR,ABS(wave)^2,time,period, $
;       XSTYLE=1,XTITLE='Time',YTITLE='Period',TITLE='Noise Wavelet', $
;       YRANGE=[MAX(period),MIN(period)], $   ;*** Large-->Small period
;       /YTYPE, $                             ;*** make y-axis logarithmic
;       NLEVELS=25,/FILL
;    IDL>
;    IDL> signif = REBIN(TRANSPOSE(signif),ntime,nscale)
;    IDL> CONTOUR,ABS(wave)^2/signif,time,period, $
;           /OVERPLOT,LEVEL=1.0,C_ANNOT='95%'
;    IDL> PLOTS,time,coi,NOCLIP=0   ;*** anything "below" this line is dubious
;
;
;----------------------------------------------------------------------------
; Copyright (C) 1995-2004, Christopher Torrence and Gilbert P. Compo
;
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.
; This routine is provided as is without any express or implied warranties
; whatsoever.
;
; Notice: Please acknowledge the use of the above software in any publications:
;    ``Wavelet software was provided by C. Torrence and G. Compo,
;      and is available at URL: http://paos.colorado.edu/research/wavelets/''.
;
; Reference: Torrence, C. and G. P. Compo, 1998: A Practical Guide to
;            Wavelet Analysis. <I>Bull. Amer. Meteor. Soc.</I>, 79, 61-78.
;
; Please send a copy of such publications to either C. Torrence or G. Compo:
;  Dr. Christopher Torrence               Dr. Gilbert P. Compo
;  Research Systems, Inc.                 Climate Diagnostics Center
;  4990 Pearl East Circle                 325 Broadway R/CDC1
;  Boulder, CO 80301, USA                 Boulder, CO 80305-3328, USA
;  E-mail: chris[AT]rsinc[DOT]com         E-mail: compo[AT]colorado[DOT]edu
;----------------------------------------------------------------------------
;-

FUNCTION morlet, $ ;*********************************************** MORLET
  k0,scale,k,period,coi,dofmin,Cdelta,psi0

  IF (k0 EQ -1) THEN k0 = 6d
  n = N_ELEMENTS(k)
  expnt = -(scale*k - k0)^2/2d*(k GT 0.)
  dt = 2*!PI/(n*k(1))
  norm = SQRT(2*!PI*scale/dt)*(!PI^(-0.25))   ; total energy=N   [Eqn(7)]
  morlet = norm*EXP(expnt > (-100d))
  morlet = morlet*(expnt GT -100)  ; avoid underflow errors
  morlet = morlet*(k GT 0.)  ; Heaviside step function (Morlet is complex)
  fourier_factor = (4*!PI)/(k0 + SQRT(2+k0^2)) ; Scale-->Fourier [Sec.3h]
  period = scale*fourier_factor
  coi = fourier_factor/SQRT(2)   ; Cone-of-influence [Sec.3g]
  dofmin = 2   ; Degrees of freedom with no smoothing
  Cdelta = -1
  IF (k0 EQ 6) THEN Cdelta = 0.776 ; reconstruction factor
  psi0 = !PI^(-0.25)
; PRINT,scale,n,SQRT(TOTAL(ABS(morlet)^2,/DOUBLE))
  RETURN,morlet
END

;****************************************************************** WAVELET
FUNCTION wavelet_xu,y1,dt, $   ;*** required inputs
  S0=s0,DJ=dj,J=j, $   ;*** optional inputs
  PAD=pad,MOTHER=mother,PARAM=param, $
  VERBOSE=verbose,NO_WAVE=no_wave,RECON=recon, $
  LAG1=lag1,SIGLVL=siglvl,DOF=dof,GLOBAL=global, $   ;*** optional inputs
  SCALE=scale,PERIOD=period,YPAD=ypad, xu_tceq7=xu_tceq7,$  ;*** optional outputs
  DO_DAUGHTER=do_daughter,DAUGHTER=daughter,COI=coi, $
  SIGNIF=signif,FFT_THEOR=fft_theor, $
  OCT=oct,VOICE=voice   ;*** defunct inputs
  
  ON_ERROR,2
  r = CHECK_MATH(0,1)
  n = N_ELEMENTS(y1)
  n1 = n
  base2 = FIX(ALOG(n)/ALOG(2) + 0.4999)   ; power of 2 nearest to N

;....check keywords & optional inputs
  IF (N_ELEMENTS(s0) LT 1) THEN s0 = 2.0*dt
  IF (N_ELEMENTS(voice) EQ 1) THEN dj = 1./voice
  IF (N_ELEMENTS(dj) LT 1) THEN dj = 1./8
  IF (N_ELEMENTS(oct) EQ 1) THEN J = FLOAT(oct)/dj
  IF (N_ELEMENTS(J) LT 1) THEN J=FIX((ALOG(FLOAT(n)*dt/s0)/ALOG(2))/dj)  ;[Eqn(10)]
  IF (N_ELEMENTS(mother) LT 1) THEN mother = 'MORLET'
  IF (N_ELEMENTS(param) LT 1) THEN param = -1
  IF (N_ELEMENTS(siglvl) LT 1) THEN siglvl = 0.95
  IF (N_ELEMENTS(lag1) LT 1) THEN lag1 = 0.0
  lag1 = lag1(0)
  verbose = KEYWORD_SET(verbose)
  do_daughter = KEYWORD_SET(do_daughter)
  do_wave = NOT KEYWORD_SET(no_wave)
  recon = KEYWORD_SET(recon)
  IF KEYWORD_SET(global) THEN MESSAGE, $
    'Please use WAVE_SIGNIF for global significance tests'
    
;....construct time series to analyze, pad if necessary
  if KEYWORD_SET(xu_tceq7) then pad=0
  ypad = y1 - TOTAL(y1)/n    ; remove mean
  IF KEYWORD_SET(pad) THEN BEGIN   ; pad with extra zeroes, up to power of 2
    ypad = [ypad,FLTARR(2L^(base2 + 1) - n)]
    n = N_ELEMENTS(ypad)
  ENDIF

;....construct SCALE array & empty PERIOD & WAVE arrays
  na = J + 1                  ; # of scales
  scale = DINDGEN(na)*dj      ; array of j-values
  scale = 2d0^(scale)*s0      ; array of scales  2^j   [Eqn(9)]
  period = FLTARR(na,/NOZERO) ; empty period array (filled in below)
  wave = COMPLEXARR(n,na,/NOZERO)  ; empty wavelet array
  IF (do_daughter) THEN daughter = wave   ; empty daughter array

;....construct wavenumber array used in transform [Eqn(5)]
  k = (DINDGEN(n/2) + 1)*(2*!PI)/(DOUBLE(n)*dt)
  k = [0d,k,-REVERSE(k(0:(n-1)/2 - 1))]

;....compute FFT of the (padded) time series
  yfft = FFT(ypad,-1,/DOUBLE)  ; [Eqn(3)]

  IF (verbose) THEN BEGIN  ;verbose
    PRINT
    PRINT,mother
    PRINT,'#points=',n1,'   s0=',s0,'   dj=',dj,'   J=',FIX(J)
    IF (n1 NE n) THEN PRINT,'(padded with ',n-n1,' zeroes)'
    PRINT,['j','scale','period','variance','mathflag'], $
      FORMAT='(/,A3,3A11,A10)'
  ENDIF  ;verbose
  IF (N_ELEMENTS(fft_theor) EQ n) THEN fft_theor_k = fft_theor ELSE $
    fft_theor_k = (1-lag1^2)/(1-2*lag1*COS(k*dt)+lag1^2)  ; [Eqn(16)]
  fft_theor = FLTARR(na)
  
;....loop thru each SCALE
  FOR a1=0,na-1 DO BEGIN  ;scale
    psi_fft=CALL_FUNCTION(mother, $
      param,scale(a1),k,period1,coi,dofmin,Cdelta,psi0)
    ;print,TOTAL((ABS(psi_fft)^2))
    if KEYWORD_SET(xu_tceq7) then psi_fft=psi_fft*sqrt(n/TOTAL((ABS(psi_fft)^2)));add by Shuang Xu on 2022-05-07
    IF (do_wave) THEN $
      wave(*,a1) = FFT(yfft*psi_fft,1,/DOUBLE)  ;wavelet transform[Eqn(4)]
    period(a1) = period1   ; save period
    ;print,'TOTAL((ABS(psi_fft)^2))/n: ',TOTAL((ABS(psi_fft)^2))/n
    fft_theor(a1) = TOTAL((ABS(psi_fft)^2)*fft_theor_k)/n
    ;print,'fft_theor(a1): ',fft_theor(a1)
    IF (do_daughter) THEN $
      daughter(*,a1) = FFT(psi_fft,1,/DOUBLE)   ; save daughter
    ;print,total(abs(fft(daughter(*,a1),-1))^2)
    IF (verbose) THEN PRINT,a1,scale(a1),period(a1), $
        TOTAL(ABS(wave(*,a1))^2),CHECK_MATH(0), $
        FORMAT='(I3,3F11.3,I6)'
  ENDFOR  ;scale

  coi = coi*[FINDGEN((n1+1)/2),REVERSE(FINDGEN(n1/2))]*dt   ; COI [Sec.3g]
  
  IF (do_daughter) THEN $   ; shift so DAUGHTERs are in middle of array
    daughter = [daughter(n-n1/2:*,*),daughter(0:n1/2-1,*)]

;....significance levels [Sec.4]
  sdev = (MOMENT(y1))(1)
  fft_theor = sdev*fft_theor  ; include time-series variance
  dof = dofmin
  signif = fft_theor*CHISQR_CVF(1. - siglvl,dof)/dof   ; [Eqn(18)]
  ;signif_mean = fft_theor*CHISQR_CVF(0.5,dof)/dof 

  IF (recon) THEN BEGIN  ; Reconstruction [Eqn(11)]
    IF (Cdelta EQ -1) THEN BEGIN
      y1 = -1
      MESSAGE,/INFO, $
        'Cdelta undefined, cannot reconstruct with this wavelet'
    ENDIF ELSE BEGIN
      y1=dj*SQRT(dt)/(Cdelta*psi0)*(FLOAT(wave) # (1./SQRT(scale)))
      y1 = y1[0:n1-1]
    ENDELSE
  ENDIF
  
  RETURN,wave(0:n1-1,*)    ; get rid of padding before returning

END