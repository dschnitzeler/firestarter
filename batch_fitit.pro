;   Run Batch_fitit to apply fs.pro > Fitit to an ensemble of Monte Carlo simulations.
;   See below, 'How to set up a Monte Carlo simulation'.
;
;   INPUT
;   test: specify which test you want to apply fs.pro > Fitit to.
;         For the list of tests that have been hard-coded, see fs.pro > Lookup_test_param.
;
;   ii_start: skip the first ii_start Monte Carlo runs.
;             Default = 0.
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
;   Batch_fitit calls 'Fitit' in fs.pro. In addition to the .dat and
;   .sav file that Fitit creates for each Monte Carlo simulation, Batch_fitit
;   also creates a file 'fs.testN.mpfitfun-status.dat' that
;   contains the exit status of the MPFIT routine (if that exit status
;   is not zero). A non-zero exit status implies that something might be 
;   wrong with the fit, see the explanation at the beginning of mpfit.pro.
;   Also problems with the matrix nversion in the function fs.pro > 
;   Find_mle_q_u_ref are written to this file.
;
;
;   EXAMPLE
;   GDL> batch_fitit, test=101
;
;   HOW TO SET UP A MONTE CARLO SIMULATION
;   First, specify the properties of the simulations you want to run
;   in fs.pro > Lookup_test_param.
;   Here you can specify the frequency coverage of the simulated
;   observations, the properties of the injected model, which
;   models to fit (cmp_arr), the maximum number of model components
;   (n_cmp_max), the number of Monte Carlo realisations in the
;   ensemble (n_mc_runs) etc.
;
;   I hard-coded options for simulating a single point source, 
;   Burn slab, Gaussian |L(RM)| distribution. These are models 1, 2, 
;   and 3 in Schnitzeler (2018), and I named tests after them. 
;   For example, test 102 (the '0' in the name of a test is used as a
;   separator) is the second test where a single point source was
;   injected in the data. In test 301 a Burn slab (model 3) was
;   injected. All tests 10* are referred to as series 1, 20* are
;   series 2, and 30* are series 3.
;
;   Test 11 injects a source that emits at two RMs. The nomenclature
;   here is more complicated, because the name of each test encodes
;   the polarized flux density of the fainter of the two sources, its
;   RM, and its intrinsic polarization angle. I made it easier for you
;   to run these test: in that case, specify 'test=11', together with the
;   values for the keywords l2_test11, chi0_2_test11, and rm2_rayleigh_test11.
;   Then the program works out the correct name for the test.
;   If required, update the arrays in Lookup_series11_test_values that
;   keep track of which simulations for test 11 have been run.
;   For example, I have used
;     Lookup_series11_test_values, tested_l2_values=tested_l2_values, $
;        tested_chi0_2_values=tested_chi0_2_values, $
;        tested_rm2_rayleigh_values=tested_rm2_rayleigh_values
;
;      for aa=0,n_elements(tested_rm2_rayleigh_values)-1 do $
;        Batch_fitit,test=11,l2_test11=1000,chi0_2_test11=0, $
;          rm2_rayleigh_test11=tested_rm2_rayleigh_values[aa]
;
;   to call Batch_fitit when I ran simulations of test series 11.
;
;   You can change the range in RM that is used to find the peak
;   polarized flux density in the residual RM spectrum in 
;   fs.pro > Fill_rm_alpha_grid.
;
;   You can change the properties of the noise that is added to the
;   source models in fs.pro > Generate_mock_data.
;  
;   Now everything is set up! run batch_fitit.pro > Batch_fitit.
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

@./fs

PRO BATCH_FITIT, test=test, ii_start=ii_start, $
    l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11

    if n_elements(test) eq 0 then begin
      print,' Batch_fitit: please specify which test you want to run.'
      print,''
      print,'  for example: batch_fitit, test=101' 
      print,''
      stop
    endif
    if n_elements(ii_start) gt 0 then ii_start=long(fix(ii_start)) $
    else ii_start=0l  ; define ii_start also in the default case.

    Lookup_test_param, test, test_str=test_str, mc_sim_path=mc_sim_path, n_mc_runs=n_mc_runs, $
      l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11  ; function in fs.pro.
    if ~file_test(mc_sim_path,/directory) then spawn,'mkdir '+mc_sim_path
;   Open a file for writing the exit status of MPFIT to (if this exit
;   status is not zero). Such an exit status can indicate a problem
;   when fitting the data, see the information at the beginning of
;   batch_fitit.pro and mpfit.pro.
    openw,99,mc_sim_path+'/fs.'+test_str+'.mpfitfun-status.dat'
;
    seed=0l
    for ii=0l,n_mc_runs-1 do begin
      printf,99,'Index '+roundoff(ii)
      Fitit, test=test, mc_index=ii, seed=seed, /keep99open, skip_this_fit= (ii lt ii_start)
    endfor
;
    close,99

    suffix= (ii_start ne 0) ? ' with ii_start='+roundoff(ii_start) : ''
    print,'Completed processing test '+roundoff(test)+suffix
END
