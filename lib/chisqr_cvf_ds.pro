;   For the chi squared distribution, find the parameter x that gives
;   you a probability P(X <= x) = prob. When you call this program,
;   you need to specify both 'prob' and the number of degrees of freedom
;   of the chi squared distribution, 'n_dof'.
;
;   Use that for the chi squared distribution, CDF(a,b) = IGAMMA(b/2,a/2).
;
;   INPUT:
;   prob : probability that you want to calculate for the chi squared
;          distribution
;
;   n_dof : the number of degrees of freedom
;
;   Note: if 'prob' is the probability you are interested in, then you
;   should call the function described in this program as follows:
;     result = Chisqr_cvf_ds(prob, n_dof)
;   The original (IDL) implementation should be called as follows:
;     result = Chisqr_cvf(1-prob, n_dof)
;
;   OUTPUT
;   Chisqr_cvf_ds returns the value of 'x' for which P(X <= x) = prob,
;   where 'prob' is specified when calling the function.
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


FUNCTION CHISQR_CVF_DS, prob, n_dof
    if prob lt 0 or prob gt 1 then $
      Message,'Chisqr_cvf_ds: the specified probability must lie between 0 and 1.'
 
    x_min=  0d  ; starting point
    x_max= 10d  ; guess
;   Increase x_max if required:
    while Igamma(n_dof/2d , x_max/2d) lt prob do x_max*=2

;   At this point, you can assume that the correct value for x lies
;   between x_min and x_max. Now use bisection to 'home in' on the 
;   correct value for x.    
    eps = 1d-10  ; precision with which you want to approximate prob
    x_new = (x_min+x_max)/2d
    val = Igamma(n_dof/2d , x_new/2d)
    while abs(val - prob) gt eps do begin
      if val lt prob then x_min=x_new else x_max=x_new
      x_new = (x_min+x_max)/2d
      val = Igamma(n_dof/2d , x_new/2d)
    endwhile

    return, x_new
END
