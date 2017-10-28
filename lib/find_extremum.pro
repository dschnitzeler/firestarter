;   Determine the location of a maximum or minimum of a data set by
;   fitting a parabola to the (x,y) data points. Use the worked-out
;   equations that minimise the vertical distance^2 between the y
;   data values and the fit.
;   mpfit often converged to a wrong solution (?), and in the case of
;   a parabola it's possible to work out the coefficients that
;   find the least-squares solution.
;
;   NOTE: all errors in y_i are assumed to be equal!
;
;   INPUT: 2 arrays, x_arr and y_arr
;   OUTPUT: the coefficients a,b and c that describe the parabola fit
;           of the form f(x) = a + c*(x-b)^2:
;
;
;   METHOD:
;   If the parabola that you're fitting has the form f(x) = a + bx + cx^2,
;   then the least-squares distance between the data and the fit is 
;   minimal when SUM(y-f(x))^2 = SUM(y-(a + bx + cx^2))^2 is minimal.
;
;   Minimising this quantity gives the following 3 equations:
;   dSUM/da = SUM 2*(y - (a + bx + cx^2))*-1 = 0, so 
;             SUM (y - (a + bx + cx^2)) = 0  (1)
;   dSUM/db = SUM (y - (a + bx + cx^2))*x = 0  (2)
;   dSUM/dc = SUM (y - (a + bx + cx^2))*x^2 = 0  (3)
;   
;   Eqn. (1) can be re-written as
;     SUM y_i      = a* SUM(1)     + b * SUM(x_i)   + c * SUM(x_i^2)
;   and Eqn. (2) and (3) as
;     SUM x_i*y_i  = a* SUM(x_i)   + b * SUM(x_i^2) + c * SUM(x_i^3)
;     SUM x_i^2*y_i = a* SUM(x_i^2) + b * SUM(x_i^3) + c * SUM(x_i^4)
;   : these are 3 equations with 3 unknowns (since you know all the
;   x_i and y_i, and you can therefore calculate the above sums)
;    
;   In matrix form:
;   (SUM(1)      SUM(x_i)    SUM(x_i^2))     (a)     (  SUM(y_i)    )
;   (SUM(x_i)    SUM(x_i^2)  SUM(x_i^3))  *  (b)  =  ( SUM(x_i*y_i) )
;   (SUM(x_i^2)  SUM(x_i^3)  SUM(x_i^4))     (c)     (SUM(x_i^2*y_i))
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


FUNCTION FIT_LEAST_SQ_SOLUTION, $
    x_arr, y_arr, silent=silent, return_all_coefficients=return_all_coefficients, nointerrupt=nointerrupt

;   To improve the numerical stability of the algorithm centre the
;   grid on av_x to reduce having to deal with large x:
    n_pts=n_elements(x_arr)
    av_x=total(x_arr,/double)/n_pts
    x_arr-=av_x
;   To further improve the stability of the algorithm (avoid
;   calculating factors of x^{2,3,4} if x is large) rescale x_arr
;   by dividing by factors of 10:
    scale_x=10^(double(floor(alog10(max(abs(x_arr))))))
    x_arr/=scale_x
    if ~keyword_set(silent) then print,' rescaling x: ',scale_x

    x1_tot = total(x_arr,/double)
    x2_tot = total(x_arr^2,/double)
    x3_tot = total(x_arr^3,/double)
    x4_tot = total(x_arr^4,/double)

    matrix=[[n_pts, x1_tot,x2_tot], $
            [x1_tot,x2_tot,x3_tot], $ 
            [x2_tot,x3_tot,x4_tot]]

    data_vector=[total(y_arr,/double), $
                 total(x_arr*y_arr,/double), $
                 total(x_arr^2*y_arr,/double)]

;   Us 'machar(/double)' to select when matrix elements are smaller than machine precision:
    eps_machine=Machar(/double) 
    eps_machine=10^(double(ceil(alog10(max(abs([eps_machine.eps,eps_machine.epsneg]))))))
    sel=where(abs(matrix)/n_pts lt eps_machine, n_sel) 
    if n_sel gt 0 then matrix[sel]=0d
;   Reduce the numbers in the matrix before inverting:
;   Use matrix^-1 = (matrix/N_pts)^-1/N_pts
    inv_matrix=Invert(matrix/n_pts,status,/double)/n_pts
    if status ne 0 then begin
      print,' Warning: matrix inversion unreliable..' 
      print,' Matrix:'
      print, matrix/n_pts
      print,' Inverse matrix:'
      print, inv_matrix
      if keyword_set(nointerrupt) then $
        return,[!VALUES.F_NAN,!VALUES.F_NAN,!VALUES.F_NAN] $
      else stop
    endif
    sel=where(abs(inv_matrix) lt eps_machine, n_sel) 
    if n_sel gt 0 then inv_matrix[sel]=0d

    goto,skip_tests
    print,'11------------'
    print,matrix/n_pts  ; this is the actual matrix that is being inverted
    print,'22------------'
    print,inv_matrix*n_pts  ; this is the actual result from the matrix inversion
;    skip_tests:
    if ~keyword_set(silent) then begin
      print,'33------------'
      print,matrix##inv_matrix  ; should be equal to the identity matrix
      print,'44------------'
    endif
    skip_tests:
    solution=inv_matrix##data_vector  
    solution=reform(solution[0,*]) 
    solution=Check_precision(solution,eps_machine)

;   Check the accuracy of the inversion, by multiplying the inverse
;   matrix with the original matrix. With infinite accuracy this
;   should produce the identity matrix.
    goto,skip_test
    tmp_arr = data_vector # matrix
    index_arr=indgen(3)
    tmp_arr(index_arr,index_arr)=0
    if max(abs(tmp_arr)) gt 1e-5 or total(abs(tmp_arr),/double)/6d gt 1e-5 then begin
      print,'  fit_extremum > fit_least_sq_solution: matrix inversion produces'
      print,'  large off-diagonal elements:'
      print,''
      print,tmp_arr
      print,''
      print,'  stopping...'
      stop
    endif
    skip_test:

;   The fit solution is for the function f(x') = a' + b'*x' + c'*x'^2,
;   where x' = (x-av_x)/scale_x. Work back to f(x) = a + b*x + c*x^2:
    aa= solution[0] - solution[1]*av_x/scale_x + solution[2]*(av_x/scale_x)^2
    bb= solution[1]/scale_x - 2*(solution[2]/scale_x)*(av_x/scale_x)
    cc= solution[2]/scale_x^2
;   Re-write the result as f(x) = a" + c"*(x-b")^2
    coeff_arr=solution ;holds a copy of the fitted parameters
    coeff_arr[0]= aa-0.25d*bb^2/cc
    coeff_arr[1]= -0.5d*bb/cc
    coeff_arr[2]= cc

    return, keyword_set(return_all_coefficients) ? coeff_arr : [coeff_arr[1], coeff_arr[0]]
;   in the latter case: return the x-coordinate of the extremum, and its value
END
