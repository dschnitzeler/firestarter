;   Add spaces to the beginning of a specified string to make it into
;   the desired length; I thought this is useful to have when making
;   tables with right-aligned columns, in particular in combination
;   with `roundoff' to make strings out of numbers (use roundoff(..,0)).
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


FUNCTION FILL_STRING, string, desired_length
    super_string='                                                   '

    n_elt=n_elements(string)
    if n_elt eq 1 then $
      return,strmid(super_string,0,desired_length-strlen(string))+string $
    else begin
      arr=strarr(n_elt)
      for ii=0l, n_elt-1 do $
        arr[ii]=strmid(super_string,0,desired_length-strlen(string[ii]))+string[ii]
      return,arr
    endelse  
END

FUNCTION FILL_STR, string, desired_length
    return,FILL_STRING(string, desired_length)
END
