;   Determine the spectral index and offset of a power law that is fit
;   to Stokes I flux densities as a function of frequency.
;   The program assumes that the frequency channels are statistically
;   independent, and that the noise in each channel is described by a
;   Gaussian distribution. The noise variances do not have to be the
;   same in all channels.
;
;   Requires MPFIT and MPFITFUN, written by Craig Markwardt. You can
;   download these from
;
;     http://cow.physics.wisc.edu/~craigm/idl/fitting.html
;
;   INPUT
;   x_arr: array containing the observed frequencies
;
;   y_arr: array containing the measured flux densities
;
;   y_err: (1-sigma) measurement uncertainties in the values stored in y_arr
;
;   x_ref: reference frequency for the power law. By default, the median of
;          x_arr is calculated using the IDL function median( ,/even).
;
;   start_point: a 2D array containing an initial guess for the flux
;          density that the source emits at the reference frequency
;          (first element) and the coefficient of the power law
;          (second element). If 'start_point' is not specified, the
;          program assumes a coefficient of -0.7 and calculates the
;          maximum likelihood estimator for the flux density at the
;          reference frequency using an analytical expression.
;
;   sp_index: set this keyword if you are fitting a power law to a
;             spectral energy distribution. This limits the range of
;             allowed spectral indices to [-6,3], which covers the
;             spectral indices of all known pulsars and most (all?)
;             active galactic nuclei (see Lorimer et al. 1995 and
;             Bates et al. 2013).
;
;   OUTPUT
;   The function 'Fit_powerlaw' returns a 2D array, which stores the
;   two parameters that describe the power law (flux density at the
;   reference frequency and the coefficient of the power law). These
;   parameters are derived by minimizing the chi squared/maximizing
;   the likelihood. The 2D array that is returned by 'Fit_powerlaw'
;   has the same structure as the array that can be specified with the
;   keyword 'start_point' when calling the function 'Fit_powerlaw'.
;
;   EXAMPLE
;   fit_result=Fit_powerlaw(freq_arr,stokes_i_arr,noise_i_arr)
;
;
;   HISTORY
;   16.11.2017 (DS) Updated Fit_powerlaw so that also spectral indices
;                   in the power-law fit to Stokes I are limited to
;                   the range [-6,3].                 
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


FUNCTION POWERLAW, x, p
    return, p[0]*x^p[1]
END

FUNCTION FIT_POWERLAW, x_arr, y_arr, y_err, x_ref=x_ref, start_point=start_point, sp_index=sp_index
    if n_elements(x_ref) eq 0 then x_ref=Median(x_arr,/even,/double)  ; the reference frequency for the power law
    if n_elements(start_point) eq 0 then begin
;     If no starting point for the Levenberg-Marquardt algorithm has
;     been defined by the user, create one:
      start_point=dblarr(2)
      start_point[1]=-0.7d  ; typical value for the power law in extragalactic radio sources
;     For this spectral index, find the maximum likelihood estimator
;     for the multiplicative (scale) factor of the power law, which is
;     called 'p[0]' in the function 'Powerlaw' above.
      start_point[0]= $
        total(y_arr*(x_arr/x_ref)^start_point[1]/y_err^2,/double)/$
        total((x_arr/x_ref)^(2*start_point[1])/y_err^2,/double)
    endif

    par_properties= replicate({fixed:0,limited:[0,0],limits:[0D,0D]}, 2)
    if keyword_set(sp_index) then begin
;     In this case a power law is fitted to a spectral energy
;     distribution. Limit the range of allowed values for the 
;     spectral index to [-6,3]:
      par_properties[1].limited[0:1]=1
      par_properties[1].limits=[-6D,3D]
    endif

    x_ratio_arr=(x_arr/x_ref)
    result=Mpfitfun('Powerlaw', x_ratio_arr, y_arr, y_err, start_point, parinfo=par_properties, $
                     covar=covar, perror=perror, status=status, errmsg=errmsg, /silent)
    if status le 0 then Message,errmsg
    if status ne 1 then print,'Fit_powerlaw status: '+roundoff(status)

;   Calculate the reduced chi squared value for this fit, assuming Gaussian
;   errors with variance y_err^2:
    n_dof=n_elements(x_arr)-2  ; the number of degrees of freedom (two-parameter model)
    chi2_red= total((y_arr-Powerlaw(x_ratio_arr,result))^2/y_err^2,/double)/n_dof

    return,[result[0], perror[0], result[1], perror[1], chi2_red]
END
