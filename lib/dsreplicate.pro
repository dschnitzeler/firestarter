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


FUNCTION DSREPLICATE1, arr, n_times
;   Clone 'arr' (assumed to be a double-precision array) n_times

    n_elt=n_elements(arr)
    new_arr=dblarr(n_elt*n_times)
    for ii=0,n_elt-1 do new_arr[ii+n_elt*indgen(n_times)]=replicate(arr[ii],n_times)
    return, new_arr
END

FUNCTION DSREPLICATE2, arr, n_multiples
;   Replace each element in 'arr' (assumed to be double-precision)
;   n_multiples times.

    n_elt=n_elements(arr)
    new_arr=dblarr(n_elt*n_multiples)
    for ii=0,n_elt-1 do $ 
      new_arr[ii*n_multiples : (ii+1)*n_multiples-1]=replicate(arr[ii],n_multiples)
    return, new_arr
END
