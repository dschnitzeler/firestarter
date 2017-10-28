;   Given the frequency coverage of the observations, calculate the
;   FWHM of the RMSF, the FWHM of any extended RM feature that you are
;   sensitive to, the maximum RM for which your observations are
;   sensitive, and RM_Rayleigh (DS 2018). All RMs are expressed in
;   units of rad m^{-2}.
;
;   INPUT:
;   freq_start, freq_end, channel_width (all in MHz)
;
;   OUTPUT:
;   FWHM RMSF, FWHM extended RM structure, largest RM, RM_Rayleigh (all in rad/m2)
;
;
;   MIT License
;
;   Copyright (c) 2017 Dominic H.F.M. Schnitzeler
;
;   Permission is hereby granted, free of charge, to any person obtaining a copy
;   of this software and associated documentation files (the "Software"), to deal
;   in the Software without restriction, including without limitation the rights
;   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
;   copies of the Software, and to permit persons to whom the Software is
;   furnished to do so, subject to the following conditions:
;
;   The above copyright notice and this permission notice shall be included in all
;   copies or substantial portions of the Software.
;
;   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
;   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
;   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
;   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
;   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
;   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
;   SOFTWARE.


FUNCTION RMSF_CHARS, freq_start, freq_end, channel_width, verbose=verbose, accurate=accurate
    if n_params() eq 0 then begin
      print,' The syntax for RMSF_chars: '
      print,' rmsf_chars(freq_start, freq_end, channel_width, verbose=verbose, accurate=accurate'
      print,' Frequencies are in MHz, RMs in rad m^(-2).'
      stop
    endif

    if freq_end lt freq_start then swap,freq_start,freq_end
    if ~keyword_set(verbose) then verbose=0

    delta_lambda2= (299.792458d/freq_start)^2-(299.792458d/freq_end)^2
    av_freq=0.5d*(freq_start+freq_end)
    small_delta_lambda2 = 2*299.792458d^2*channel_width/av_freq^3*(1+0.5d*(channel_width/av_freq)^2)  
;   = equation 35 in Brentjens & De Bruyn (2005)

    fwhm_rmsf=3.8d/delta_lambda2
    widest_rm=!dpi/(299.792458d/freq_end)^2
;   Calculate the largest RM that you can measure according to Brentjens & De Bruyn:
    largest_rm=1.9d/small_delta_lambda2 
;   Note that this is somewhat wider than the sqrt(3)/small_delta_lambda2 from B&dB (2005)
;   The referee for Schnitzeler et al. (2009) was very helpful, and
;   provided us with this equation, which he explained was more
;   accurate than the equation from Brentjens & De Bruyn (2005).
;
;   RM_Rayleigh is a crude way for estimating if a source is resolved,
;   similar to the Rayleigh criterion in optics (see equation 14 in 
;   Schnitzeler 2018):
    rm_rayleigh = !dpi/((299.792458d/freq_start)^2-(299.792458d/freq_end)^2)
    
;   Check that RM_Rayleigh is small enough that the frequency channels
;   with top-hat response functions in frequency can be approximated
;   as channels with top-hat response functions in wavelength squared:
;   (= Equation 14 from Schnitzeler & Lee 2015)
    rm_crit = 14400d*(freq_start/1d3)^2.5d*(freq_end/1d3)^0.5d/channel_width

    if keyword_set(verbose) then begin
      print,'Characteristics of the RMSF:'   
      print,'  FWHM of the RMSF: '+roundoff(fwhm_rmsf,1)+' rad/m2'
      print,'  FWHM of widest detectable RM structure: '+roundoff(widest_rm,1)+' rad/m2'
      print,'  largest detectable RM: '+roundoff(largest_rm,1)+' rad/m2'
      print,'  boundary RM (top-hat response): '+roundoff(rm_crit,1)+' rad/m2'
      if rm_rayleigh gt rm_crit then suffix='*' else suffix=''
      print,'  RM_Rayleigh: '+roundoff(rm_rayleigh,1)+suffix+' rad/m2'
    endif

    return,[fwhm_rmsf,widest_rm,largest_rm,rm_rayleigh]
END
