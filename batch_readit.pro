;   Run Readit on an entire ensemble of Monte Carlo simulations, and
;   write the results to the file called 'fs.testN.summary.dat'.
;
;   INPUT
;   test: the identification for the test, as defined in 
;         fs.pro > Lookup_test_param
;
;   If you want to run Batch_readit on a double point source model,
;   use test=11 and specify the following parameters to describe the
;   second source:
;   l2_test11: polarized flux density of the second source component
;
;   chi0_2_test11: intrinsic polarization angle \chi_0 of the second
;                  source component.
;
;   rm2_rayleigh_test11: value for the RM of the second source
;                        component, in units of RM_Rayleigh.
;   No defaults.
;
;
;   OUTPUT
;   A file 'fs.testN.summary.dat' is written to the directory
;   that stores the data for the Monte Carlo simulations. Each line in
;   this file shows the ranking of models according to the Akaike
;   Information Criterion, the Bayesian Information Criterion, and the
;   reduced chi squared, and model-averaged parameter values as follows:
;
;    ranking according to AIC (AIC weight of each fit) |
;
;    model-averaged parameters, using the AIC weights  ||
;
;    ranking according to BIC (BIC weight of each fit) |
;
;    model-averaged parameters, using the BIC weights  ||
;
;    ranking according to the reduced chi squared (reduced chi
;     squared for each fit) 
;
;   The ID that is listed for a model refers to the row number+1 of that
;   model in the list that is created by fs.pro > fill_model_type_grid.
;   For example, ID = "1" is the first model in this list etc.
;   The number between round brackets is the weight that is calculated
;   from the value of the AIC or BIC metric (see Schnitzeler 2018 and
;   references therein).
;
;   The order in which the model-weighted parameters are listed is the
;   same as the order in which parameters are written to the file
;   'fs.dat' by fs.pro > Fitit, i.e.,
;   (1) RM (= the RM_0 or RM_c from section 2.1 in Schnitzeler 2018)
;   (2) sigma_RM_intern
;   (3) delta_RM
;   (4) sigma_RM_extern
;   (5) the spectral index \alpha of the polarized flux density
;       spectrum of the source. Note that the program uses a power law
;       to describe this intrinsic emission of the source (intrinsic =
;       before the emission is depolarized). The spectral index
;       \alpha, together with the flux densities in Stokes Q and U
;       at the reference frequency (columns 6 and 7) define this power law.
;   (6) Stokes Q at the reference frequency, Q_ref
;   (7) Stokes U at the reference frequency, U_ref.
;       Q_ref and U_ref, together with \alpha, describe the intrinsic
;       polarized flux density spectrum of the source.
;   If n_cmp_max is larger than 1, then all these parameters are
;   specified for each source component that was fitted, without
;   using separators like "|".
;
;
;   EXAMPLE
;   GDL> batch_readit, test=101
;
;
;   HOW TO SET UP A MONTE CARLO SIMULATION
;   See batch_fitit.pro.
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

@./readit

PRO BATCH_READIT, $
    test=test, l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11

    if n_elements(test) eq 0 then begin
      print,' Batch_readit: please specify the test you want to run.'
      print,''
      print,'  for example: batch_readit, test=101' 
      print,''
      stop
    endif
    write_summary_file = 1  ; boolean, write information to the file fs.testN.summary.dat (value = 1) or not.
    silent = 1  ; boolean, write information to terminal (silent = 0) or not.

    Lookup_test_param, test, test_str=test_str, mc_sim_path=mc_sim_path, n_mc_runs=n_mc_runs, $
      l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11  ; function in fs.pro
    mc_files=File_search(mc_sim_path+'/fs.'+test_str+'.index*.dat', count=n_mc_runs_tmp)
    if n_mc_runs ne n_mc_runs_tmp then begin
      print,' test '+roundoff(test)
      print,' Batch_readit: n_mc_runs is not equal to the number of files from the Monte Carlo simulation.'
      print,' '+roundoff(n_mc_runs_tmp)+' files were found, there should be '+roundoff(n_mc_runs)+' files.'
      print,' stopping..'
      stop
    endif

    print,' Reading data for test '+roundoff(test)

    for ii=0l,n_mc_runs-1 do begin
      if ii/25. eq fix(ii/25.) then print,ii
      Readit, test=test, mc_index=ii, silent=silent, summary=summary  ; function in readit.pro

      if keyword_set(write_summary_file) then begin   
        if ii eq 0 then openw,11,mc_sim_path+'/fs.'+test_str+'.summary.dat'  ; Initialize   
        printf,11,summary
      endif
    endfor

    if keyword_set(write_summary_file) then close,11
END
