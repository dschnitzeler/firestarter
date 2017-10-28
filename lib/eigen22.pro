;   Calculate the eigenvalues and eigenvectors for a 2D real matrix,
;   using the fast method explained on
;     http://www.math.harvard.edu/archive/21b_fall_04/exhibits/2dmatrices/index.html
;     (downloaded: '~/notes/Eigenvalues and eigenvectors of 2x2 matrices')
;   which is explained on
;     https://math.stackexchange.com/questions/395698/fast-way-to-calculate-eigen-of-2x2-matrix-using-a-formula
;     (downloaded: '~/notes/linear algebra - Fast way to calculate Eigen of a 2x2 matrix using a formula')
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


PRO EIGEN22, matrix, val1=val1, val2=val2, vec1=vec1, vec2=vec2, eps_machine=eps_machine
    if n_elements(eps_machine) eq 0 then eps_machine=Find_eps_machine()

;   First calculate the two eigenvalues:
    tt= matrix[0,0] + matrix [1,1]  ; the trace of the matrix
    det= matrix[0,0]*matrix[1,1] - matrix[1,0]*matrix[0,1]  ; the determinant of the matrix
;
    arg_tmp = tt^2/4d - det
;   Due to limited numerical precision, arg_tmp can get close to
;   zero. Check if this is the case:
    if arg_tmp gt eps_machine then begin
      val1= tt/2d + sqrt(arg_tmp)  ; the first eigenvalue
      val2= tt/2d - sqrt(arg_tmp)  ; the second eigenvalue
    endif $
    else begin
      val1= tt/2d
      val2= tt/2d
    endelse

;   Now calculate the eigenvectors:
    if abs(matrix[0,1]) gt eps_machine then begin
      vec1= [val1-matrix[1,1], matrix[0,1]]
      vec2= [val2-matrix[1,1], matrix[0,1]]
      return
    endif
    if abs(matrix[1,0]) gt eps_machine then begin
      vec1= [matrix[1,0], val1-matrix[0,0]]
      vec2= [matrix[1,0], val2-matrix[0,0]]
      return
    endif
;   In this case both matrix[1,0] and matrix[0,1] are zero, which
;   means that the matrix is a diagonal matrix. Then the eigenvectors
;   are simply:
    vec1= [1,0]
    vec2= [0,1]    
END

