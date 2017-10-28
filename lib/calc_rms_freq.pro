;   Calculate the RM spectrum for an array of Stokes Q and U as a
;   function of frequency, using the formalism from Schnitzeler & Lee (2015)
;   or Brentjens & De Bruyn (2005).
;
;   INPUT
;   Stokes_q: array of measured Stokes Q flux densities as a function of frequency
;             Can by a single line of sight, or a data cube.
;
;   Stokes_u: array of measured Stokes U flux densities as a function of frequency
;             Can by a single line of sight, or a data cube.
;
;   Freq_arr: frequencies at which Stokes Q and U were measured (MHz)
;             Note: the program can handle irregularly sampled data.
;
;   Freq_step: channel width (MHz)
;
;   RM_min: smallest RM for which you want to calculate the RM spectrum (rad m^-2)
;
;   RM_max: largest RM for which you want to calculate the RM spectrum (rad m^-2)
;
;   RM_step: interval in RM for calculating the RM spectrum (rad m^-2)
;
;
;   OUTPUT
;   vec_P_RM: contains the (complex) RM spectrum for the specified
;             range in RMs. Real_pat(vec_P_RM) will give you the values
;             of Stokes Q as a function of RM, Imaginary(vec_P_RM)
;             the values of Stokes U.
;
;   OPTIONS (KEYWORDS)
;   use_new_formalism: use the formalism from DS & Lee (2015) instead
;                      of Brentjens & De Bruyn (2005). This is the
;                      default. It requires calculating the error
;                      function of a complex number, which requires
;                      IDL (in Oct. 2017 it had not been implemented
;                      yet in GDL).
;   use_old_formalism: force the program to use the formalism by
;                      Brentjens & De Bruyn (2005).
;   level:             'level' times RM_critical specifies when you want to   
;                      switch from B05 to the new formalism. If it isn't set  
;                      it is assumed to be equal to 1. RM_critical is
;                      defined by equation 14 in DS & Lee (2015).
;   return_freq_spectrum: for a source that emits at a single RM, return
;                         Stokes Q and U as a function of frequency.
;                         Currently only implemented for use in
;                         combination with DS & Lee (2015), not for
;                         Brentjens & De Bruyn (2005).
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


FUNCTION INT_NEG, freq, rm
;   Returns the complex number that is calculated by equation 13 in DS
;   & Lee (2015): this is the number you need to derotate the
;   polarization vectors that have been measured as a function of frequency.
;
;   Works for both scalar and vector arrays 'freq'

    lambda = 299.792458d/freq
    phi = 2*rm*lambda^2
    return, dcomplex(cos(-phi),sin(-phi))*freq + $
      sqrt(abs(rm)*!dpi)*299.792458d*dcomplex(sign(rm),1)*erf(sqrt(abs(rm))*lambda*dcomplex(sign(rm),1))
END


FUNCTION INT_POS, freq, rm
;   Use this function to predict a frequency spectrum for a source that emits at a single rm.
    return, Int_Neg(freq, -rm)
END


FUNCTION CALC_RMS_FREQ,$
    stokes_q, stokes_u, freq_arr, freq_step, RM_min, RM_max, RM_step, $
    use_new_formalism=use_new_formalism,use_old_formalism=use_old_formalism, $
    level=level, return_freq_spectrum=return_freq_spectrum

    if n_elements(level) eq 0 then level=1d

;   Set up the (central) wavelength squared array:
    lambda2_arr=0.5d*(299.792458d/freq_arr)^2 + 0.5d*(299.792458d/(freq_arr+freq_step))^2
;   see DS & Lee (2015) for this definition of lambda^2_central
;   Note: you need freq_step, not abs(freq_step) here because the mean
;   lambda2_arr depends on whether the frequency increases or
;   decreases across the channel.

    n_RM=ceil(1d*(RM_max - RM_min)/RM_step)+1 
;   Note: RM_min + n_RM*RM_step must be >= RM_max
    RM_arr=dindgen(n_RM)*RM_step+RM_min

    n_freq= n_elements(freq_arr)
    freq_min=min(freq_arr,max=freq_max)
;   Evaluate equation 14 from DS & Lee (2015) to determine if you can
;   use the formalism by Brentjens & De Bruyn (2005), or need to use 
;   our new formalism:
    RM_crit = 14400d*(freq_min/1d3)^2.5d*(freq_max/1d3)^0.5d/abs(freq_step)  ; in rad/m2

;   Save computing time when running our new RM synthesis algorithm:
;   Calculate the integral Int_Neg (equation 12 from the DS & Lee 2015)
;   for the low-end frequencies of all channels. Then check if you can
;   use already calculated values of Int_Neg for the high-end
;   frequencies of the channels.
;
;   Assume that frequencies have already been ordered in ascending or
;   descending order. Check which of the two:
    delta_freq_asc = freq_arr[1]-freq_arr[0] gt 0
    freq_arr_high= freq_arr+abs(freq_step)  ; the high-end frequency of each channel
    if delta_freq_asc then $
      sel=where(abs(shift(freq_arr,-1) - freq_arr_high) gt abs(freq_step)/100d, n_sel, compl=sel_compl, ncompl=n_sel_compl) $
    else $
      sel=where(abs(shift(freq_arr, 1) - freq_arr_high) gt abs(freq_step)/100d, n_sel, compl=sel_compl, ncompl=n_sel_compl) 
;   'sel' contains indices for which you have to calculate Int_Neg,
;   while for indices in 'sel_compl' you can re-use values that you
;   have calculated already.
;   (Another way of seing this: indices in 'sel' refer to
;   non-contiguous channels, while those in 'sel_compl' refer to
;   contiguous channels.)

;   Treat the way that you calculate the RM spectrum of one line of
;   sight in the same way as how you calculate the RM spectrum of a
;   data cube. Therefore, if you're dealing with a single line
;   of sight, cast it into the same 3D format as that of a data cube:
    dimen=size(stokes_q,/dim)
    n_dimen=n_elements(dimen)
    case n_dimen of 
    1: begin
         stokes_q=reform(stokes_q,1,1,dimen[0])
         stokes_u=reform(stokes_u,1,1,dimen[0])

         dimen=[1,1,dimen[0]]
       end
    2: begin
         stokes_q=reform(stokes_q,dimen[0],1,dimen[1])
         stokes_u=reform(stokes_u,dimen[0],1,dimen[1])

         dimen=[dimen[0],1,dimen[1]]
       end
    else: ctu=1
;     in that case you're dealing with Stokes Q/Stokes U data cubes
    endcase

;   Construct a vec_P(freq) and a vec_P(RM) data cube:
    tmp=dblarr(dimen[0],dimen[1],n_RM)
    vec_p=dcomplex(stokes_q,stokes_u) 
    vec_p_RM=dcomplex(tmp,tmp)
    lambda2_cube=rebin(reform(lambda2_arr,1,1,n_freq),dimen[0],dimen[1],n_freq)

;   The actual calculation of the RM spectrum. 
    for ii=0l, n_RM-1 do begin
      RM=RM_arr[ii]  ; calculate the RM spectrum at this RM

      if (abs(RM) lt level*RM_crit or keyword_set(use_old_formalism)) and ~keyword_set(use_new_formalism) then begin
;       Use Brentjens & De Bruyn (2005)
        if keyword_set(return_freq_spectrum) then begin
          print,'return_freq_spectrum is not implemented yet in calc_rms_freq...' & stop
        endif

        arg=-2*RM*lambda2_cube 
        vec_p_RM[*,*,ii]=1d/n_freq*total(vec_p*dcomplex(cos(arg),sin(arg)),3,/double)
      endif $
      else begin
;       Use our RM synthesis formalism and frequency channels with a
;       top-hat channel response in frequency.

        if !version.os eq 'darwin' then begin
;         'erf' cannot handle complex arguments in GDL, therefore
;         write out an error:
          print,' Calc_rms_freq: erf cannot handle complex arguments in GDL.'
          print,' Please use IDL!' & close,/all & stop
        endif

;       Step 1: calculate the net derotation vectors for the low-end frequencies of all frequency channels:
        int_freq_one = Int_Neg(freq_arr, RM)
;       Step 2: look up or calculate the high-end frequencies of all frequency channels:
        int_freq_two = int_freq_one*0
        if delta_freq_asc gt 0 then begin
          if n_sel gt 0 then int_freq_two[sel]= Int_Neg(freq_arr[sel]+freq_step, RM)
;         Note: you need freq_step, not abs(freq_step) here because
;         the 'mean derotation vector' depends also on nu_2 = nu_1 +
;         freq_step, not nu_1 + abs(freq_step).
          if n_sel_compl gt 0 then int_freq_two[sel_compl]= int_freq_one[sel_compl+1]
;           The high-end frequency of the current channel is equal to max(nu_1,nu_2) of the next channel.
        endif $  ; if delta_freq_asc gt 0
        else begin
          if n_sel gt 0 then int_freq_two[sel]= Int_Neg(freq_arr[sel]+freq_step, RM)
          if n_sel_compl gt 0 then int_freq_two[sel_compl]= int_freq_one[sel_compl-1]
;           The high-end frequency of the current channel is equal to max(nu_1,nu_2) of the previous channel
        endelse
        p_derot = (int_freq_two-int_freq_one)/freq_step  ; the channel-average derotation 
;       freq_step, not abs(freq_step): the net derotation vector for a
;       single channel should not depend on whether nu_1 {<,>} nu_2 (the
;       frequencies of the channel edges). If nu_1 > nu_2, freq_step =
;       nu_2 - nu_1 is negative, and cancels the minus sign of the
;       integral which you get when you change the direction of integration.
;       Therefore you get the same result as when nu_1 < nu_2.
        if keyword_set(return_freq_spectrum) then begin
;         For a source that emits at a single RM, return Stokes Q and U as a function of frequency:
          return,transpose([[freq_arr],[Real_Part(p_derot)],[Imaginary(p_derot)]])
        endif
        p_derot /= (abs(p_derot)*n_freq)  ; turn into unit vectors to conserve flux density, and normalise (='/n_freq')
        vec_p_RM[*,*,ii] = total(vec_p*p_derot,/double)
;
      endelse
    endfor  ; for ii=0l, n_RM-1

;   If you're calculating the RM spectrum of 1 line of sight,
;   or only 1 row of lines of sight, then create the right output 
;   format for vec_p_RM:
    if n_dimen eq 1 or n_dimen eq 2 then vec_p_RM=reform(vec_p_RM)

;   Return the (complex) RM synthesis spectrum, one line per RM:
    return,vec_p_RM
END
