;   Convert 'number' to a string with 'no_digits' digits behind the
;   decimal point.
;
;   EXAMPLE
;   print,roundoff(3.14159,2)
;   output: 3.14
;
;   NOTE
;   Roundoff works up to and including 18 significant figures.
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

FUNCTION ROUNDOFF,number,no_digits, use_d=use_d
;   Use /use_d to replace the decimal dot with a 'd'

    n_numbers=n_elements(number)

    if n_params() eq 1 then $
;     Assume that only 'number' has been specified:     
      no_digits=0

    tmp_number_arr=strarr(n_numbers)

    for ii=0l,n_numbers-1 do begin
      if no_digits gt 0 and number[ii] ne 0d then begin
;       when a number is exactly 0 then the alog10(number) will
;       produce a meaningless number. Therefore, separate numbers == 0,
;       and make 'roundoff' write them out as 0, which is done by
;       treating '0' in the same way as when no_digits=0

        tmp_number=double(number[ii])*10d^no_digits
        tmp_number=floor(tmp_number+0.5d,/l64)/10d^no_digits  ; the actual formatting
        if tmp_number eq 0 then begin
          format='(f'+strcompress(2+no_digits,/rem)+'.'+strcompress(no_digits,/rem)+')'
          tmp_number=strcompress(string(tmp_number,format=format),/rem)
          goto,skip_this_zero
        endif
;       required when converting tmp_number to string format:
        tmp=floor(alog10(abs(tmp_number)))  ; the number of digits before the decimal '.'
        tot_n_digits=no_digits+1+ tmp*(tmp gt 0)+1
;       the number of digits in the number that is requested is the
;       number of digits behind the decimal '.' + 1 digit for the
;       decimal '.' itself + the number of digits before the decimal
;       dot.
;       The last '+1' is required because every number will have at 
;       least one digit in front of the decimal '.', when alog10 is
;       still < 1, and tmp therefore still 0.
;       add another digit if the original number is smaller than 0:
        tot_n_digits+=(sign(number[ii]) lt 0)

;       strcompress(tmp_number,/rem) doesn't work, you really
;       need to explicitly format the string first, and then use 
;       strcompress:
        format='(f'+strcompress(tot_n_digits,/rem)+'.'+strcompress(no_digits,/rem)+')'
;       '+1' since the decimal '.' in tmp_number also takes up one character
        tmp_number=strcompress(string(tmp_number,format=format),/rem)

        skip_this_zero:
        tmp_number_arr[ii]=tmp_number   
      endif $  ; if no_digits gt 0 and number[ii] ne 0d
      else $
        tmp_number_arr[ii]=strcompress(floor(number[ii]+0.5d,/l64),/rem)

      if keyword_set(use_d) then begin
;       Replace the decimal dot with a 'd'
        tmp=Strsplit(tmp_number_arr[ii],'.',/extract)
        tmp2= tmp[0]+'d'+tmp[1]
        tmp_number_arr[ii]=tmp2
      endif
    endfor ;for i=0,n_numbers-1

    if n_numbers eq 1 then return,tmp_number_arr[0] $
    else return,tmp_number_arr
END

