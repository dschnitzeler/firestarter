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


FUNCTION CHECK_PRECISION, var, eps_machine
    if n_params() eq 1 then begin
;     Check if 'var' is float or double:
      eps_machine=Machar(double= (Size(var,/type) eq 5)) 
      eps_machine=10^(double(ceil(alog10(max(abs([eps_machine.eps,eps_machine.epsneg]))))))
    endif

    if n_elements(var) eq 1 then begin
      if abs(var) lt eps_machine then var=0d
    endif $
    else begin
      sel=where(abs(var) lt eps_machine, n_sel)
      if n_sel gt 0 then var[sel]=0d
;     assuming 'var' is double precision to start with
    endelse
    return,var
END
