;   Firestarter is a QU fitting-based algorithm for extracting
;   polarization information from noisy data sets. The algorithm is
;   described in Schnitzeler (2018). Please refer to this paper when
;   you use the program. There are two main differences compared to
;   competing algorithms available at the time of writing:
;   (1) the spectral indices of polarized sources are fitted, instead
;       of assuming they are scaled versions of Stokes I, followed by
;       dividing Stokes Q and U by Stokes I (assuming the sources emit 
;       a power law spectrum, typical for synchrotron emission, before 
;       this emission gets depolarized), and
;   (2) the program includes the properties of the noise in Stokes Q,
;       U, and I when fitting source models to the data. The noise in 
;       these Stokes parameters is assumed to follow Gaussian distributions
;       with zero mean. The standard deviations of these distributions 
;       in Q and U do not have to be the same, and are allowed to vary 
;       across the frequency band.
;   Frequency channels, and measurements of Stokes Q and U in each
;   channel, are assumed to be statistically independent.
;
;   To apply the Firestarter algorithm to real data, or to mock
;   observations from a Monte Carlo simulation, run fs.pro > Fitit.
;   For example, if you want to analyse your data, you can run
;
;     Fitit, meas_arr=meas_arr, cmp_arr=[1,2], n_cmp_max=2
;
;   which will get you pretty far: you tell the main fitting program,
;   'Fitit', which data to use, the types of models being fitted and
;   the maximum number of source components. These variables are explained 
;   below, in the section 'INPUT'. Alternatively, to find out what
;   these three parameters mean, you can simply type 'Fitit' at the
;   command prompt in your terminal, without specifying anything else.
;   The program will then write a short description of each parameter
;   and the layout of the array meas_arr to the terminal. 
;
;   For batch processing (for example, if you want to run Fitit on a
;   large number of Monte Carlo simulations) use batch_fitit.pro > 
;   Batch_fitit. To make the output from Fitit easier to interpret, 
;   use readit.pro > Readit and, if required, batch_readit.pro >
;   Batch_readit (see the explanation on what these programs do at the
;   start of each program).
;
;   The program is written in the Gnu Data Language (GDL), an
;   open-source version of the Interactive Data Language (IDL). It is
;   compatible with GDL version 0.9.6 and IDL version 8.1.
;
;   Comments / bugs reports (or bug fixes) / feedback are welcome! 
;   E-mail me at d<insert my surname>(at)<name of that e-mail service
;   provided by Google>.com.
;
;   The algorithm was developed for analysing data from the 16 cm band
;   on the Australia Telescope Compact Array, which covers frequencies
;   between 1.3-3.1 GHz (nominally down to 1.1 GHz, but these low
;   frequencies are normally not usable because of radio-frequency
;   interference). If you want to run your own Monte Carlo
;   simulations, have a look at 'how to set up a Monte Carlo
;   simulation' in batch_fitit.pro to see which variables can be adjusted.
;
;   The significance of the detection is calculated using the log
;   likelihood ratio, the signal-to-noise ratio by comparing the
;   extent of the 1-sigma error ellipse in the (Q_ref,U_ref) plane
;   with the polarized flux density of the signal you detected.
;   See appendix A from Schnitzeler (2018). 
;   The program calculates the signal-to-noise ratio only for the
;   highest peak in the residual RM spectrum (the RM spectrum that is
;   calculated from the measured Stokes Q and U flux densities minus
;   the best-fitting model prediction for these flux densities).
;   As pointed out in section 6.3 from Schnitzeler (2018), if the
;   noise variances in Q and U vary with frequency, then the peak
;   with the highest polarized flux density is not automatically the 
;   peak with the highest signal-to-noise ratio. One could calculate 
;   the signal-to-noise ratio of the polarization vector at each RM, 
;   but this slows down the program.
;
;   You can speed up processing by changing the parameters 'rm_min',
;   'rm_max', and 'rm_step' in the function Fill_rm_alpha_grid below. 
;   Changing rm_min and rm_max also allows you to search for sources
;   with very large positive or negative RMs.
;
;   If you want to run a Monte Carlo simulation, you can change the
;   noise variances in Stokes Q and U in 'Generate_mock_data' below.
;
;   You might be wondering if Firestarter is a fancy acronym. It
;   isn't.. By sharing this software I hope you will have more time to 
;   do your own science, and make the most of your data. If you're
;   using this code for that purpose, or if you want to build a better 
;   algorithm, either way, I hope you become more enthusiastic about 
;   radio polarimetry. Hence, 'Firestarter'. 
;
;   A few basic rules that will help you read the program code:
;   v It's easy to jump between these procedures by searching for the
;     comment string ";'".
;   v Variable names ending in "_tmp" have a very local scope.
;   v Nomenclature:
;     'n_cmp_max' specifies the maximum number of model components
;         that should be used in a fit.
;         Models with up to n_cmp_max components are fitted to the data.
;     'p' is the array used by MPFITFUN, which only contains
;         non-linear model parameters.
;     'p_big' = 'p' + the linear model parameters Q_ref and U_ref
;     'p_arr' is only used in Simple_Gaussian and Simple_IFD, and is
;         a subarray of p_big with a length of n_param_max cells, which
;         describes a single source component.
;     'fit_result' is the best fitting model, which is produced by MPFITFUN.
;         This array includes the linear model parameters.
;   v Tests follow the naming convention <cmp_list>0<index>:
;     <cmp_list> indicates which model components are injected. 
;       One can choose between the eight models listed in section 2.1
;       in Schnitzeler (2018).
;     '0' is a separator, which makes it easier for the user to find the ID of the test,
;     'index' specifies which test was run using cmp_list.
;     For example: 
;      test 201 injects one or more Gaussian sources (= model 2) to the data. 
;      test 301 uses a Burn slab instead of a Gaussian source (= model 3).
;   v DS18 = Schnitzeler (2018).
;   v ML = maximum likelihood
;     LM = Levenberg-Marquardt (algorithm)
;     logl = log likelihood
;     loglr = log likelihood ratio
;     chi2_red = reduced chi squared
;     chi2_red_max = maximum reduced chi squared for a model fit,
;                    see section 4 in DS18.
;     AIC = Akaike Information Criterion
;     BIC = Bayesian Information Criterion
;     det_significance = significance of the detection
;
;   And, finally: CHECK YOUR RESULTS to see if they make sense!
;
;
;   INPUT
;   meas_arr: array with seven columns, and a number of rows equal to
;     the number of frequency channels. The seven columns are:
;     (0) frequency (in MHz) 
;     (1,2) Stokes I flux density and its 1-sigma uncertainty (both in mJy)
;     (3,4) Stokes Q flux density and its 1-sigma uncertainty (both in mJy)
;     (5,6) Stokes U flux density and its 1-sigma uncertainty (both in mJy)
;     No defaults.
;
;    cmp_arr: lists which model components can be fitted. Eight different
;      components are defined, see figure 2 and section 2.1 in DS18.
;      If Firestarter is applied to real data, the default choice is
;      model type 1, a point source in RM. In a Monte Carlo simulation, all
;      eight model types are included (see the function Lookup_test_param).
;
;    n_cmp_max: models with up to and including n_cmp_max components
;      are fitted to the data. If Firestarter is applied to real data,
;      by default n_cmp_max = 5. For Monte Carlo simulations, the value 
;      depends on which test was run (see the function Lookup_test_param).
;
;    src: name of the source being analysed (not required).
;
;    silent: boolean. If set to 1, no information is written to the
;      terminal and no data are plotted. By default = 0.
;
;    keep99open: boolean. By default, a new file
;      'fs.mpfitfun-status.dat' (see 'OUTPUT' below) is created each
;      time the program is run. By using '/keep99open', the file
;      stream is not closed when Fitit has finished, which means that
;      consecutive runs of Fitit all write to the same file
;      'fs.mpfitfun-status.dat'. For example, batch_fitit.pro >
;      Batch_fitit opens file stream 99 once, then calls fs.pro >
;      Fitit repeatedly using the keyword '/keep99open'. 
;
;   Use the following keyword to create mock observations as part of a
;   Monte Carlo simulation:
;    test: specifies which test to analyse (see the list of tests in
;          Lookup_test_param below). Default is test = 100.
;
;    mc_index: specifies which Monte Carlo run from 'test' to analyse.
;              mc_index starts at zero, default = 0.
;
;    seed: seed used in the random number generator RANDOMN. default = 0l.
;
;    skip_this_fit: boolean, default = 0. Used in combination with
;      batch_fitit.pro > batch_fitit when ii_start is not equal to zero:
;      this skips the first ii_start fits.
;
;   If you want to run Fitit, but not Batch_fitit, on a double point
;   source model, then use test=11 and specify the following parameters:
;    l2_test11: polarized flux density of the second source component.
;
;    chi0_2_test11: intrinsic polarization angle \chi_0 of the second
;                   source component.
;
;    rm2_rayleigh_test11: value for the RM of the second source
;                         component, in units of RM_Rayleigh.
;   No defaults.
;
;
;   OUTPUT
;   fs.pro > Fitit writes information on the model fits to the terminal
;   and to a number of files on your local hard drive. It also plots
;   the data and the model fits in a window. 
;
;   fs.pro > Fitit reports information on the model fits to the
;   terminal window, and it writes this information to a file that is
;   called 'fs.<src name>.dat' or 'fs.dat', depending on
;   whether you have specified a value for the keyword 'src'.
;   If you run a Monte Carlo simulation, the file is called
;   'fs.testN.indexN.dat', where 'N' is an integer number with one or
;   more digits.
;   Each line that is written to the terminal or to the file describes
;   the results from fitting a single source model, using this layout:
;
;   <which model was fitted> | <metrics> | <source parameters and their uncertainties>
;
;   If models with more than one component were fitted, then a
;   vertical bar is used to distinguish the set of parameters of each
;   of the source components. Instead of
;
;   <source parameters and their uncertainties>
;
;   the file contains 
;
;   <parameters of component 1 and their uncertainties> |
;   <parameters of component 2 and their uncertainties> | (etc.).
;   
;   The numbers in front of the first "|" indicate which model was
;   fitted to the data, and refer to the numbers listed in section 2.1 
;   from DS18 and in the function Lookup_model below. There can be
;   more than one number if you're fitting models with more than one
;   component (n_cmp_max > 1 in that case).
;   E.g., "1 1" indicates two sources that each emit at a single RM.
;
;   The numbers between the first and the second "|" are the metrics
;   of the fit. These are:
;   (1) the log likelihood ratio
;   (2) the Akaike Information Criterion (AIC)
;   (3) the Bayesian Information Criterion (BIC)
;   (4) the reduced chi squared of the fit
;   (5) the maximum value for the reduced chi squared of the fit.
;       If number (4) is larger than number (5), the differences
;       between the model fit and the data are so large that they are
;       unlikely to be the result of the noise: this indicate a poor
;       fit to the data. See also the discussion in the penultimate
;       paragraph in section 4 of DS18.
;   (6) the significance of the detection, calculated from the log
;       likelihood ratio, and translated to the equivalent detection
;       significance for a 1D Gaussian distribution, to make the number
;       easier to interpret. The number "5" is listed even if
;       the source has been detected at a higher significance.
;       This has to do with the way the detection significance is
;       calculated from the log likelihood ratio, and numerical
;       accuracy.
;   (7) flag status. The flag status can take on the following values:
;       flag = 0: everything ok
;       flag = 1,2,3,4: problem with finding the uncertainties in
;                       Q_ref or U_ref in Find_err_q_u_ref.
;       flag = 100: problem with finding the uncertainties in Q_ref or
;                   U_ref in Find_err_q_u_ref.
;       flag = 111: highest peak in the residual RM spectrum has a 
;                   signal-to-noise ratio that is below the threshold 
;                   (= you’re fitting noise). 
;       flag = 1000: problem with matrix inversion in Find_mle_q_u_ref. 
;       flag = 2000: MPFIT reached the maximum number of iterations
;                    (MPFIT status 5). This can indicate that MPFIT
;                    did not converge on the right parameter solution.
;                    Note: the status of MPFIT is written to the file 
;                    fs.mpfitfun-status.dat or
;                    fs.testN.mpfitfun-status.dat if MPFIT exits with
;                    a status other than 1. An exit status <= 0 forces
;                    Fitit to stop automatically.
;       Any flag status other than 0 indicates a problem with fitting
;       the data; don't trust these fits. Therefore I discard them
;       from further analysis (e.g., ranking models based on their
;       AIC, BIC, etc.)
;
;   The numbers after the second "|" list the parameters of the
;   best-fitting model together with their uncertainties. 
;   The order of these parameters is always the same:
;   (1,2) RM and its uncertainty
;   (3,4) sigma_RM_intern and its uncertainty (see section 2.1 in DS18).
;   (5,6) delta_RM and its uncertainty 
;   (7,8) sigma_RM_extern and its uncertainty 
;   (9,10) the spectral index \alpha that was fitted to the
;          polarization data, and its uncertainty. The program assumes
;          that the source emits intrinsically (before any
;          depolarization) a power law. This power law is
;          characterized by a spectral index \alpha, and by the flux
;          densities in Stokes Q and U at the reference frequency
;          (columns 11-16).
;   (11-13) the value for Stokes Q at the reference frequency, Q_ref,
;           that, together with the spectral index (parameter 9)
;           describe the intrinsic emission of the source in Stokes Q,
;           before this emission is depolarized.
;           Note that two uncertainties are listed: one is the
;           uncertainty extending towards smaller (more negative)
;           values of Q_ref, the other is the uncertainty towards larger
;           values of Q_ref. This comes from the fitting of a parabola
;           to determine the uncertainty in Q_ref (and U_ref).
;   (14-16) the value for Stokes U at the reference frequency, U_ref,
;           before this emission is depolarized (the equivalent of the
;           numbers listed in 11-13, but for Stokes U instead of
;           Stokes Q).
;   Note that parameters that were not fitted to the data (e.g.,
;   delta_RM when fitting an external Gaussian depolarizing screen)
;   are assigned a value of zero; the same is true of the
;   uncertainties of these parameters.
;
;   readit.pro > Readit helps you interpret these numbers.
;
;   Alternatively, if you apply fs.pro > Fitit to a Monte Carlo
;   simulation, then Fitit creates two files for each mock data set:
;   (1) 'fs.testN.indexN.dat' ('N' is an integer).
;       This file has the same layout as 'fs.<src name>.dat' or 'fs.dat'.
;   (2) 'fs_mock_data.testN.indexN.sav'
;       This file stores the data (Stokes Q and U as a function of
;       frequency) that were generated for each ensemble of the Monte
;       Carlo simulation. It is an IDL sav file, which can be opened
;       using the command 'restore' (see the IDL documentation).
;   For example, when analysing the Monte Carlo simulations for
;   DS18 the program created files with names like
;   'fs.test101.index1.dat' and
;   'fs_mock_data.test101.index1.sav'
;
;   Fitit also plots the data and the model fits. More information
;   about this can be found in the function Plotit below.
;
;   Fitit creates a file that is named either fs.mpfitfun-status.dat
;   or fs.testN.mpfitfun-status.dat (with N the number of the test,
;   see Lookup_test_param below). These files show the exit status of
;   MPFIT (if this exit status is not equal to zero), which can help
;   identify problems when fitting the data. An exit status <= 0
;   forces Fitit to stop automatically. Also problems with the matrix
;   inversion in the function Find_mle_q_u_ref are written to these 
;   files.
;   If a problem occurs, then a short description of the model being
;   fitted is written to the file. The numbers used there refer to the 
;   list of model in section 2.1 of DS18. 
;   For example, '1 1 1 1 1 status: 2' means that MPFIT exited with 
;   status = 2 when fitting a model consisting of five point sources
;   in RM.
;
;
;   EXAMPLE
;    GDL > Fitit, meas_arr=meas_arr, src='0823-500', cmp_arr=[1,2], n_cmp_max=2, /silent
;
;
;   HOW TO SET UP A MONTE CARLO SIMULATION
;   See batch_fitit.pro.
;
;   HISTORY
;   06.11.2017 (DS) Updated Plotit to fix a problem when plotting coloured, filled circles in IDL
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

@./lib/fill_string
@./lib/find_extremum
@./lib/dsreplicate
@./lib/rmsf_chars
@./lib/roundoff
@./lib/find_eps_machine
@./lib/calc_rms_freq
@./lib/inv22
@./lib/eigen22
@./lib/check_precision
@./lib/sign
@./lib/echo
@./lib/find_plot_range
@./lib/chisqr_cvf_ds

PRO LOOKUP_MODEL, $
    model_type, type_str=type_str, param_used=param_used, n_param_used=n_param_used, $
    param_names_full=param_names_full, n_param_max=n_param_max

;   Given a certain model type, find the name of the model, the name
;   of the function used to predict Stokes Q and U as a function of
;   frequency, and whether a variable is a free parameter or is kept fixed.
;   Parameters are stored in the following order:
;    (0) RM, (1) sigma_RM_intern, (2) \Delta_RM, (3) sigma_RM_extern,
;    (4) \alpha, (5) Q_ref, (6) U_ref

    if n_params() eq 0 then begin
      param_names_full = ['RM','sigma_RM_intern','delta_RM','sigma_RM_extern','alpha','Q_ref','U_ref']
      n_param_max = n_elements(param_names_full)  ; the maximum number of free parameters a model can have.
      return
    endif

    Case model_type of 
      1: begin
           type_str='Point source'
           param_used=[1,0,0,0,1,1,1]
         end
      2: begin
           type_str='Gaussian external to the source'
           param_used=[1,0,0,1,1,1,1]
         end
      3: begin
           type_str='Burn slab'
           param_used=[1,0,1,0,1,1,1]
         end
      4: begin
           type_str='Burn slab + Gaussian external to the source'
           param_used=[1,0,1,1,1,1,1]
         end
      5: begin
           type_str='Internal Faraday dispersion'
           param_used=[1,1,1,0,1,1,1]
         end
      6: begin
           type_str='Internal Faraday dispersion + Gaussian external to the source'
           param_used=[1,1,1,1,1,1,1]
         end
      7: begin
           type_str='Gaussian turbulence inside the source'
           param_used=[1,1,0,0,1,1,1]
         end
      8: begin
           type_str='Gaussian turbulence inside source + Gaussian external to the source'
           param_used=[1,1,0,1,1,1,1]
         end
    endcase

    n_param_used=total(param_used)
END

PRO FILL_FNC_NAME_ARR, model_type=model_type, fnc_name=fnc_name
;   Fill the array 'fnc_name' (created before you call 'Fill_fnc_name_arr').

    sel=where(model_type ge 0 and model_type le 2,n_sel,compl=sel_compl,ncompl=n_sel_compl)
    if n_sel ge 1 then fnc_name[sel]='Simple_Gaussian'
    if n_sel_compl ge 1 then fnc_name[sel_compl]='Simple_IFD'
END

FUNCTION FIND_MODEL_TYPE_ARR_INPUT, test
;   Find out which model or models was/were injected in the mock data

    tmp=Strsplit(Roundoff(test),'0',/extract)
    n_cmp_tmp=strlen(tmp[0])  ; the number of source components in the injected model
    model_type_arr_input=intarr(n_cmp_tmp)
    for ii=0,n_cmp_tmp-1 do model_type_arr_input[ii]=Fix(Strmid(tmp[0],ii,1))

    return, model_type_arr_input
END

PRO FILL_RM_ALPHA_GRID, rm_half, fwhm_rmsf
    common rm_alpha_grid, rm_min, rm_max, rm_step, alpha_min, alpha_max, alpha_step

;   To find starting values for the Levenberg-Marquardt algorithm, Fitit
;   calculates RM spectra between rm_min and rm_max in intervals of rm_step.
;
;   rm_half is the 'maximum RM' as defined by Brentjens & De Bruyn
;     (2005; see also the discussion in Schnitzeler & Lee 2015). 
;     rm_half is in units of rad m^{-2}.
;   rm_step is the full-width at half maximum of the RM spread
;     function (in rad m^{-2}). 
;
;   I haven't written a program for calculating RM spectra for \alpha not equal
;   to zero, therefore, keep alpha_min=alpha_max=0.

    rm_max=rm_half & rm_min=-rm_max & rm_step = fwhm_rmsf/4D 
    alpha_min=0D & alpha_max=0D & alpha_step=1D
END

PRO LOOKUP_SERIES11_TEST_VALUES, $
    tested_l2_values=tested_l2_values, tested_chi0_2_values=tested_chi0_2_values, $
    tested_rm2_rayleigh_values=tested_rm2_rayleigh_values

;   Which models did you test so far?
    tested_l2_values = [1000D, 100D, 10D]  ; polarized flux density of the second source
    tested_chi0_2_values = [0D, 45D, 90D, 135D]  ; intrinsic angle of the second source (in degrees)
    tested_rm2_rayleigh_values = [5,3,1,0.8,0.6,0.4,0.2]  ; value of RM_2, in units of RM_Rayleigh
END

PRO FIND_DOUBLE_POINT_SOURCE_MODEL_INDICES, $
    l2=l2, chi0_2=chi0_2, rm2_rayleigh=rm2_rayleigh, nr_for_test=nr_for_test, str_for_test=str_for_test, $
    translate_index=translate_index

;   Specify the polarized flux density of the second peak, and its RM
;   and reduced chi squared, then work out the index of the test where
;   this model was run.
;   The idea is that this program calculates the correct index for
;   you, instead of you having to check the index manually.
;   This program stores both the index for the test (as a 64-bit long)
;   and the test as a text string, using the variables 'test' and 'test_str'.
;
;   If you use the keyword 'translate_index' then nr_for_test or
;   str_for_test is translated into the values for L2, chi0_2, and RM2
;   for that test.
;
;   Note: rm2 is expressed in units of RM_Rayleigh.

;   Which models did you test so far?
    Lookup_series11_test_values, tested_l2_values=tested_l2_values, $
      tested_chi0_2_values=tested_chi0_2_values, $
      tested_rm2_rayleigh_values=tested_rm2_rayleigh_values

    if ~keyword_set(translate_index) then begin
;     Check if the model you're querying has been simulated:
      sel1=where(tested_l2_values eq l2,n_sel1)
      sel2=where(tested_chi0_2_values eq chi0_2,n_sel2)
      sel3=where(tested_rm2_rayleigh_values eq rm2_rayleigh,n_sel3)
      if n_sel1*n_sel2*n_sel3 eq 0 then begin
;       In this case at least one of sel1,sel2,sel3 is zero.
        print,' Find_double_point_source_model_indices: the model you specified'
        print,' has not been specified yet... Care to run it now?'
        close,/all & stop
      endif
 
      str_for_test='110'+roundoff(sel1[0]+1)+'0'+roundoff(sel2[0]+1)+'0'+roundoff(sel3[0]+1)
      nr_for_test=Long64(str_for_test)
    endif $
    else begin
      if n_elements(str_for_test) eq 0 then $
;       In this case assume that 'nr_for_test' is set:
        str_for_test=roundoff(nr_for_test)
;     Assume that either str_for_test or nr_for_test is specified.      

      ind_arr=Strsplit(Strmid(str_for_test,3),'0',/extract)
      l2= tested_l2_values[ind_arr[0]-1]
      chi0_2= tested_chi0_2_values[ind_arr[1]-1]
      rm2_rayleigh= tested_rm2_rayleigh_values[ind_arr[2]-1]
    endelse
END

PRO LOOKUP_TEST_PARAM, $
    test, test_str=test_str, mc_sim_path=mc_sim_path, cmp_arr=cmp_arr, $
    n_cmp_max=n_cmp_max, n_mc_runs=n_mc_runs, $
;   The following three variables are only important when you are
;   testing a double point source model:
    l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11

;   Define parameters for each of the tests.

    common input_model, model_type_arr_input, fnc_name_arr_input, p_big_input1, p_big_input2

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common l_injected, l1, l2

    common rm_alpha_grid, rm_min, rm_max, rm_step, alpha_min, alpha_max, alpha_step

;   --- Specify parameters relating to the observations ---
    freq_min = 1300D  ; MHz
    freq_step = 8D
    n_chan = 229l  
    freq_ref = 2100D  ; reference frequency for the power spectrum in flux density
    freq_arr=freq_min+dindgen(n_chan)*freq_step
    freq_max=freq_min+n_chan*freq_step
    lambda2_arr=(299.792458D/freq_arr)^2

;   Define parameters for sampling the (RM,\alpha) grid:
    tmp= Rmsf_chars(freq_min, freq_max, freq_step)
    fwhm_rmsf=tmp[0]
    rm_half=tmp[2]  ; the 'largest detectable RM' from Brentjens & De Bruyn (2005).
;                     see also Schnitzeler & Lee (2015).
    Fill_rm_alpha_grid, rm_half, fwhm_rmsf
    rm_rayleigh=tmp[3]

;   Specify where on your hard drive the Monte Carlo simulations
;   should be stored:
    tmp=Strsplit(test,'0',/extract)
    mc_sim_path='fs_mc_sim_series'+roundoff(tmp[0])

;   --- Specify parameters that define the injected model ---
;   Define parameters for two source components, which you then can
;   update in the Case statement. 
    model_type_arr_input=Find_model_type_arr_input(test)
    n_cmp_inj_max=n_elements(model_type_arr_input)
;
;   Generic Component 1:
    l1 = 10D  ; Flux density of source 1 at the reference frequency, arbitrary units 
    chi0_1 = 0D  ; angle in degrees
    rm1 = 0D  ; rad m^-2
    sigma_RM_intern1 = 0D & delta_RM1 = 0D & sigma_RM_extern1 = 0D  ; in rad m^-2
    alpha1 = 0D
;
;   Generic Component 2:
    l2 = 0D
    chi0_2 = 0D
    rm2 = 0D
    sigma_RM_intern2 = 0D & delta_RM2 = 0D & sigma_RM_extern2 = 0D
    alpha2 = 0D

    tmp_str=roundoff(test)
    if strmid(tmp_str,0,2) eq '11' then begin
;     Specify test=11 if you want to analyse a double point source model.
      l1 = 1000D  ; ensures that this source is much stronger than the noise
;      n_cmp_max = 1  
      n_cmp_max = 2

      if strlen(tmp_str) eq 2 then begin
;       This is the first time you call Lookup_test_param for a
;       double point source model, therefore you need to generate the
;       index for the test and convert this number to a string.
        if (n_elements(l2_test11) eq 0 or n_elements(chi0_2_test11) eq 0 or n_elements(rm2_rayleigh_test11) eq 0) then begin
          print,' Lookup_test_param: for test 110*, please specify L2, chi0_2 and/or RM2'
          print,'   using the keywords l2_test11, chi0_2_test11, and rm2_rayleigh_test11.'
          close,/all
          stop
        endif 
;       Replace the value of L2, chi0_2, and RM2 with the
;       values specified when calling Lookup_test_param:
        l2=l2_test11 & chi0_2=chi0_2_test11 & rm2_rayleigh=rm2_rayleigh_test11
        rm2=rm2_rayleigh_test11*rm_rayleigh
        Find_double_point_source_model_indices, $
          l2=l2, chi0_2=chi0_2, rm2_rayleigh=rm2_rayleigh, nr_for_test=nr_for_test, str_for_test=str_for_test
        test=nr_for_test & test_str=str_for_test
      endif $  ; if strlen(tmp_str) eq 2
      else begin
;       In this case you're testing a double point source
;       model, and already know the correct index for the test
;       (parameters nr_for_test and str_for_test).
;       Look up the parameters for the second source:
        Find_double_point_source_model_indices, /translate_index, $
          l2=l2, chi0_2=chi0_2, rm2_rayleigh=rm2_rayleigh, nr_for_test=test, str_for_test=str_for_test               
        rm2=rm2_rayleigh*rm_rayleigh
        test_str=str_for_test
      endelse

      mc_sim_path+='_n_cmp_max_'+roundoff(n_cmp_max)+'/test'+str_for_test
    endif  ; if strmid(tmp_str,0,2) eq '11'

;   ---- Specify various other parameters ---
    n_mc_runs = 10  ; how many Monte Carlo realisations do you need?
    cmp_arr  = [1,2,3,4,5,6,7,8]  ; which model components do you want to test?
cmp_arr = 1 ;;
;   The maximum number of source components you want to search is
;   specified by the parameter n_cmp_max in the 'case' statement below.
    test_str= 'test'+roundoff(test)

    if strmid(roundoff(test),0,2) eq '11' then goto,skip_case_statement
;   In the Case statement you can modify parameters for source components 1 and 2.
    Case test of
     100: begin
;         Default for testing the program
          n_mc_runs=100
;         Note: at the end of this program you make the array
;         model_type_arr_input have n_cmp_max columns (the reason for
;         this is explained in the comments later on).
          n_cmp_max = 1  ; search sources with at most this many components
          end
     101: begin
          n_cmp_max = 1  
          end
     102: begin
          n_cmp_max = 2
          end
;    Tests 103 and 104 choose as starting point for MPFITFUN sigma_RM_intern=0, \delta_RM=0, sigma_RM_extern=0
     103: begin
          n_cmp_max = 1  
          fwhm_rmsf *= -1
          end
     104: begin
          n_cmp_max = 2
          fwhm_rmsf *= -1
          end
;    Tests 105 and 106 use the same setup as tests 101/102, but reduce the signal-to-noise ratio of the source.
     105: begin
          l1 = 1D
          n_cmp_max = 1  
          end
     106: begin
          l1 = 1D
          n_cmp_max = 2  
          end
;    Reduce the signal-to-noise ratio even further:
     107: begin
          l1 = 0.1D
          n_cmp_max = 1  
          end
     108: begin
          l1 = 0.1D
          n_cmp_max = 2  
          end
;    Compare to a model with no signal at all (the control sample):
     109: begin
          l1 = 0D
          n_cmp_max = 1  
          end
    1010: begin
          l1 = 0D
          n_cmp_max = 2  
          end
;
     201: begin
          l1=50D
          sigma_RM_extern1 = rm_rayleigh
          n_cmp_max = 1  
          end
     202: begin
          l1=50D
          sigma_RM_extern1 = rm_rayleigh/4
          n_cmp_max = 1  
          end
     203: begin
          l1=50D
          sigma_RM_extern1 = rm_rayleigh/25
          n_cmp_max = 1  
          end
;    Fit for two components:
     204: begin
          l1=50D
          sigma_RM_extern1 = rm_rayleigh
          n_cmp_max = 2  
          end
     205: begin
          l1=50D
          sigma_RM_extern1 = rm_rayleigh/4
          n_cmp_max = 2  
          end
     206: begin
          l1=50D
          sigma_RM_extern1 = rm_rayleigh/25
          n_cmp_max = 2  
          end
;    Reduce the signal-to-noise ratio:
     207: begin
          l1=5D
          sigma_RM_extern1 = rm_rayleigh
          n_cmp_max = 1  
          end
     208: begin
          l1=5D
          sigma_RM_extern1 = rm_rayleigh/4
          n_cmp_max = 1  
          end
     209: begin
          l1=5D
          sigma_RM_extern1 = rm_rayleigh/25
          n_cmp_max = 1  
          end
;    Again fit also for two components:
     2010: begin
          l1=5D
          sigma_RM_extern1 = rm_rayleigh
          n_cmp_max = 2  
          end
     2011: begin
          l1=5D
          sigma_RM_extern1 = rm_rayleigh/4
          n_cmp_max = 2  
          end
     2012: begin
          l1=5D
          sigma_RM_extern1 = rm_rayleigh/25
          n_cmp_max = 2  
          end
;    Reduce the signal-to-noise ratio even further:
     2013: begin
          l1=0.5D
          sigma_RM_extern1 = rm_rayleigh
          n_cmp_max = 1  
          end
     2014: begin
          l1=0.5D
          sigma_RM_extern1 = rm_rayleigh/4
          n_cmp_max = 1  
          end
     2015: begin
          l1=0.5D
          sigma_RM_extern1 = rm_rayleigh/25
          n_cmp_max = 1  
          end
;    And again fit also for two components:
     2016: begin
          l1=0.5D
          sigma_RM_extern1 = rm_rayleigh
          n_cmp_max = 2  
          end
     2017: begin
          l1=0.5D
          sigma_RM_extern1 = rm_rayleigh/4
          n_cmp_max = 2  
          end
     2018: begin
          l1=0.5D
          sigma_RM_extern1 = rm_rayleigh/25
          n_cmp_max = 2  
          end
;
     301: begin
          l1=25D
          delta_RM1 = 2*rm_rayleigh 
          n_cmp_max = 1  
          end
     302: begin
          l1=25D
          delta_RM1 = rm_rayleigh
          n_cmp_max = 1  
          end
     303: begin
          l1=25D
          delta_RM1 = rm_rayleigh/4
          n_cmp_max = 1  
          end
     304: begin
          l1=25D
          delta_RM1 = rm_rayleigh/25
          n_cmp_max = 1  
          end
;    Fit for two components:
     305: begin
          l1=25D
          delta_RM1 = 2*rm_rayleigh
          n_cmp_max = 2  
          end
     306: begin
          l1=25D
          delta_RM1 = rm_rayleigh
          n_cmp_max = 2  
          end
     307: begin
          l1=25D
          delta_RM1 = rm_rayleigh/4
          n_cmp_max = 2  
          end
     308: begin
          l1=25D
          delta_RM1 = rm_rayleigh/25
          n_cmp_max = 2  
          end
;    Reduce L1:
     309: begin
          l1=2.5D
          delta_RM1 = 2*rm_rayleigh 
          n_cmp_max = 1  
          end
     3010: begin
          l1=2.5D
          delta_RM1 = rm_rayleigh
          n_cmp_max = 1  
          end
     3011: begin
          l1=2.5D
          delta_RM1 = rm_rayleigh/4
          n_cmp_max = 1  
          end
     3012: begin
          l1=2.5D
          delta_RM1 = rm_rayleigh/25
          n_cmp_max = 1  
          end
;    Fit for n_cmp_max = 2:
     3013: begin
          l1=2.5D
          delta_RM1 = 2*rm_rayleigh
          n_cmp_max = 2  
          end
     3014: begin
          l1=2.5D
          delta_RM1 = rm_rayleigh
          n_cmp_max = 2  
          end
     3015: begin
          l1=2.5D
          delta_RM1 = rm_rayleigh/4
          n_cmp_max = 2  
          end
     3016: begin
          l1=2.5D
          delta_RM1 = rm_rayleigh/25
          n_cmp_max = 2  
          end
;    Reduce L1 even further:
     3017: begin
          l1=0.25D
          delta_RM1 = 2*rm_rayleigh 
          n_cmp_max = 1  
          end
     3018: begin
          l1=0.25D
          delta_RM1 = rm_rayleigh
          n_cmp_max = 1  
          end
     3019: begin
          l1=0.25D
          delta_RM1 = rm_rayleigh/4
          n_cmp_max = 1  
          end
     3020: begin
          l1=0.25D
          delta_RM1 = rm_rayleigh/25
          n_cmp_max = 1  
          end
;    Fit for n_cmp_max = 2:
     3021: begin
          l1=0.25D
          delta_RM1 = 2*rm_rayleigh
          n_cmp_max = 2  
          end
     3022: begin
          l1=0.25D
          delta_RM1 = rm_rayleigh
          n_cmp_max = 2  
          end
     3023: begin
          l1=0.25D
          delta_RM1 = rm_rayleigh/4
          n_cmp_max = 2  
          end
     3024: begin
          l1=0.25D
          delta_RM1 = rm_rayleigh/25
          n_cmp_max = 2  
          end
;
     901: begin
;         Test if the program exits gracefully if the highest peak in
;         the residual RM spectrum has a lower flux density than the
;         cutoff specified by pi_cutoff in Fit_model.
          cmp_arr=[1]
          n_cmp_max = 3
          n_mc_runs = 100
          end
    else: begin
          print,' Test '+roundoff(test)+' has not been defined yet. Eager to add it?'
          close,/all & stop
          end
    endcase
    skip_case_statement:

;   --- After the Case statement come variables and arrays that are
;       the same for all tests ---

    n_elt_tmp=n_elements(model_type_arr_input)
    if n_elt_tmp lt n_cmp_max then begin
;     Ensure that the array model_type_arr_input has n_cmp_max columns,
;     the same number of columns as model_type_grid and model_type_arr.
      tmp=intarr(n_cmp_max) 
      tmp[0:n_elt_tmp-1]=model_type_arr_input
;     Cells with indices from n_elt_tmp to n_cmp_max-1 are kept at zero, indicating no component.
      model_type_arr_input=tmp
    endif

    if n_cmp_inj_max gt 1 then $
      if model_type_arr_input[1] ne 0 and l2 eq 0 then begin
;       Safeguard that you're being consistent:
        print,' Lookup_test_param: you specified more than one component for the injected model,'
        print,' but the second component has l2 = 0. Please check.'
        close,/all & stop
      endif

    fnc_name_arr_input=strarr(n_cmp_inj_max)
;   Both fnc_name_arr_input and model_type_arr_input have n_cmp_inj_max
;   cells at this point.
    Fill_fnc_name_arr, model_type=model_type_arr_input, fnc_name=fnc_name_arr_input

;   Store the properties of the injected sources:
    q_ref_1=l1*cos(2*chi0_1/!radeg) & u_ref_1=l1*sin(2*chi0_1/!radeg)
    q_ref_2=l2*cos(2*chi0_2/!radeg) & u_ref_2=l2*sin(2*chi0_2/!radeg)
    p_big_input1=[rm1, sigma_RM_intern1, delta_RM1, sigma_RM_extern1, alpha1, q_ref_1, u_ref_1]
    p_big_input2=[rm2, sigma_RM_intern2, delta_RM2, sigma_RM_extern2, alpha2, q_ref_2, u_ref_2]
END

PRO READ_MEASUREMENTS, $
    meas_arr, stokes_i_arr=stokes_i_arr, noise_i_arr=noise_i_arr, $
    cmp_arr=cmp_arr, n_cmp_max=n_cmp_max

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common rm_alpha_grid, rm_min, rm_max, rm_step, alpha_min, alpha_max, alpha_step

;   Extract information on the frequency coverage of the observations,
;   and the measured flux densities in Stokes I, Q, and U and the noise
;   level in Stokes V (sigma_V):
    freq_arr=reform(meas_arr[0,*])  ; [MHz]
    n_chan=n_elements(freq_arr)  ; the number of frequency channels
    freq_ref=Median(freq_arr,/even)  ; the reference frequency for the power-law fit.
    lambda2_arr=(299.792458D/freq_arr)^2
    stokes_i_arr=reform(meas_arr[1,*])
    noise_i_arr=reform(meas_arr[2,*])
    stokes_q_arr=reform(meas_arr[3,*])
    noise_q_arr=reform(meas_arr[4,*])
    stokes_u_arr=reform(meas_arr[5,*])
    noise_u_arr=reform(meas_arr[6,*])

;   Define parameters for sampling the (RM,\alpha) grid:
    freq_min=min(freq_arr,max=freq_max)
    freq_step=abs(min(freq_arr-shift(freq_arr,-1)))
    tmp= Rmsf_chars(freq_min, freq_max, freq_step)
    fwhm_rmsf=tmp[0]
    rm_half=tmp[2]  ; the 'largest detectable RM' from Brentjens & De Bruyn (2005).
;                     see also Schnitzeler & Lee (2015).
    Fill_rm_alpha_grid, rm_half, fwhm_rmsf 

;   Define which models can be fitted:
    if n_elements(cmp_arr) eq 0 then cmp_arr = 1  ; fit only point sources in RM.
    if n_elements(n_cmp_max) eq 0 then n_cmp_max = 5  ; fit models up to n_cmp_max components
END

FUNCTION CALC_CRUDE_Q_U_ERRORS, noise_arr, alpha_opt
    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

;   Estimate the uncertainty in Q_ref or U_ref in a crude way, using 
;   equation 21 in DS18, assuming that the noise variances in Stokes Q 
;   and U are equal in each channel, and that there are no other
;   sources in the RM spectrum.
;   Note that this expression for the measurement uncertainties in
;   Q_ref and U_ref reduces to the equation derived by Macquart et
;   al. (2012) if the noise variances in Stokes Q and U are equal and
;   constant across the frequency band, and you're investigating
;   sources with \alpha=0 or if you're calculating RM spectra
;   (where \alpha is zero by construction).
    return, 1d/total(1d/noise_arr^2 * (freq_arr/freq_ref)^(2*alpha_opt),/double) 
END

FUNCTION FIND_CHI2_RED_MAX, n_param_tot=n_param_tot
    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

;   Calculate the cut-off value for the reduced chi squared.

    if n_elements(n_param_tot) eq 0 then $
      Message,'Find_chi2_red_max: please specify the number of parameters using the keyword n_param_tot'

;   For this calculation to work, I'm assuming that chi squared
;   follows a Gaussian distribution, which requires n_chan >> 1. 
;   (see section 4 in DS18).
;   Check that this is the case:
    if n_chan lt 50 then $
      Message,' Readit: n_chan too small to use the current definition of chi2_red_max.'

;   Calculate the number of degrees of freedom:
    n_dof=2*n_chan-n_param_tot  
    return, 1 + 5*sqrt(2d/n_dof)
;   This ensures that the chances of finding a model with a chi2_red
;   >~~ chi2_red_max are the same, independent of n_chan.
;   The '5' in the expression for chi2_red_max refers to 5-sigma in a
;   normal distribution. 
END

PRO SELECT_RELIABLE_FITS, $
    metrics_arr, param_arr, n_param_max, chi2_red_max_override=chi2_red_max_override, $
    sel=sel, n_sel=n_sel, compl=sel_compl, ncompl=n_sel_compl

;   Check that the fitted \alpha is not -6 or +3, the min. or max. allowed value:
    dimen=size(param_arr,/dim)
    n_cmp=dimen[0]/n_param_max  ; how many model components were fitted?
    ind_arr=indgen(n_cmp)
    use_arr=intarr(dimen[1])+1  
    for ii=0,dimen[1]-1 do begin
      sel=where(param_arr[4+ind_arr*n_param_max,ii] eq -6 or param_arr[4+ind_arr*n_param_max,ii] eq 3, n_sel)
      if n_sel ne 0 then use_arr[ii]=0
    endfor

;   Check if at least one model fit is usable:
    if total(use_arr) eq 0 then begin
;     In this case all model fits converged on alpha=-6 or alpha=+3.
      sel=-1 & n_sel=0 & n_sel_compl=dimen[1] & sel_compl=indgen(n_sel_compl)
      return
    endif

;   Only select model fits with a reduced chi squared smaller than
;   chi2_red_max or chi2_red_max_override, and flag status equal to zero: 
    if n_elements(chi2_red_max_override) eq 0 then $
      sel=where(metrics_arr[4,*] le metrics_arr[5,*] and metrics_arr[7,*] eq 0 and use_arr eq 1, $  
                n_sel, compl=sel_compl, ncompl=n_sel_compl) $
    else $
      sel=where(metrics_arr[4,*] le chi2_red_max_override and metrics_arr[7,*] eq 0 and use_arr eq 1, $
                n_sel, compl=sel_compl, ncompl=n_sel_compl)
END

FUNCTION FIND_L_INJECTED_AND_N_CMP_MAX, test=test
    common l_injected, l1, l2

    Lookup_test_param, test, n_cmp_max=n_cmp_max

    return, [l1,l2,n_cmp_max]
END

FUNCTION FIND_SERIES_STR, test=test
    tmp=Strsplit(roundoff(test),'0',/extract) 

    return, tmp[0]
END

PRO FILL_MODEL_TYPE_GRID, $
    cmp_arr=cmp_arr, n_cmp_max=n_cmp_max, n_models=n_models, $
    model_type_grid=model_type_grid, fnc_name_grid=fnc_name_grid, family_tree_grid=family_tree_grid

;   Fill_model_type_grid depends on the following input parameters:
;   cmp_arr: specifies list of which model types to consider (1d array)
;   n_cmp_max: specifies the maximum number of components in a single model (scalar)

    n_diff_cmp=double(n_elements(cmp_arr))
    n_models= (n_diff_cmp gt 1) ? n_diff_cmp*(n_diff_cmp^n_cmp_max-1)/(n_diff_cmp-1) : n_cmp_max
    model_type_grid=dblarr(n_cmp_max,n_models)
    fnc_name_grid=strarr(n_cmp_max,n_models)
    family_tree_grid=intarr(2,n_models) 
;   stores the row index of the 'parent model' which is the model that has one 
;   component less than the current model, plus the flag status of the
;   parent model after Fit_model has run.

    index=0l 
    for ii=0,n_cmp_max-1 do begin
      length_this_ii= n_diff_cmp^(ii+1)

      for jj=0,ii do begin
        clone_arr= Dsreplicate2(cmp_arr,n_diff_cmp^(ii-jj))
        model_type_grid[jj, index : index+length_this_ii-1]= Dsreplicate1(clone_arr,n_diff_cmp^jj)
      endfor

      if ii ne 0 then begin
        tmp_arr=indgen(length_previous_ii)+index_prev_ii
        family_tree_grid[0,index : index+length_this_ii-1]= $
          Dsreplicate2(tmp_arr,n_diff_cmp) 
;       Note: each ii has n_diff_cmp as many entries as ii-1.
;       Therefore, fill family_tree_grid by cloning the row entries
;       from block ii-1 n_diff_cmp times.
      endif $
      else family_tree_grid[0,0:length_this_ii-1]=replicate(-1,n_diff_cmp)   
;     family_tree_grid[0,ii] eq -1 indicates that there is no parent model

      index_prev_ii=index & length_previous_ii=length_this_ii
      index+=length_this_ii
    endfor

    Fill_fnc_name_arr, model_type=model_type_grid, fnc_name=fnc_name_grid
END

PRO FILL_PARAM_USED_ARR, $
    model_type_grid=model_type_grid, n_param_max=n_param_max, n_cmp_max=n_cmp_max, $
    n_model=n_models_tested, param_used_arr=param_used_arr

;   Create the array 'param_used_arr' that tells you which parameters
;   were used (allowed to vary) in each model that was tested.

    param_used_arr=intarr(n_param_max*n_cmp_max,n_models_tested)  
    for ii=0,n_models_tested-1 do begin
      Lookup_model, model_type_grid[0,ii], param_used=param_used
      param_used_arr[0 : n_param_max-1,ii]=param_used
      if n_cmp_max gt 1 then begin
        for jj=1,n_cmp_max-1 do $
          if model_type_grid[jj,ii] ne 0 then $
;           Because of the way model_type_grid is constructed, you can
;           simply look up param_arr in the first rows of param_used_arr:
            param_used_arr[jj*n_param_max : (jj+1)*n_param_max-1,ii]= $
              param_used_arr[0 : n_param_max-1, model_type_grid[jj,ii]-1]
      endif  ; if n_cmp_max gt 1
    endfor  ; for ii=0,n_models_tested-1
END

FUNCTION CALC_NOISE_COVARIANCE, rm
;   Calculate the covariance matrix from noise_q_arr and noise_u_arr
;   at \alpha = 0. Use alpha=0 because you're searching the residual
;   RM spectrum, where alpha=0 by construction.

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common machine_precision, eps_machine

;   See the erratum on DS & Lee (2017) for the correct form of these equations.
    delta1= [cos(2*rm*lambda2_arr),sin(2*rm*lambda2_arr)]
    delta2= [-sin(2*rm*lambda2_arr),cos(2*rm*lambda2_arr)]
;    dd_inv=dblarr(2l*n_chan,2l*n_chan)
;    ind_arr_tmp=indgen(2l*n_chan)
;    dd_inv[ind_arr_tmp,ind_arr_tmp]=[1d/noise_q_arr^2,1d/noise_u_arr^2]
;   dd is a diagonal matrix, therefore store only the non-zero elements:
    dd_inv_diag= [1d/noise_q_arr^2,1d/noise_u_arr^2]

    matrix= dblarr(2,2)
;    matrix[0,0]= delta1##(dd_inv##delta1)
;    matrix[1,0]= delta1##(dd_inv##delta2)
;    matrix[0,1]= delta2##(dd_inv##delta1)
;    matrix[1,1]= delta2##(dd_inv##delta2)
;   Using dd_inv_diag saves time, this avoids multiplying a large
;   (2*n_chan,2*n_chan) diagonal matrix with a vector:
    matrix[0,0]= Total(delta1*(dd_inv_diag*delta1),/double)
    matrix[1,0]= Total(delta1*(dd_inv_diag*delta2),/double)
    matrix[0,1]= Total(delta2*(dd_inv_diag*delta1),/double)
    matrix[1,1]= Total(delta2*(dd_inv_diag*delta2),/double)

    sel=where(abs(matrix) lt eps_machine, n_sel)
    if n_sel gt 0 then matrix[sel]=0
 
    return,Inv22(matrix)  ; the covariance matrix
END

FUNCTION CALC_SNR, rm, polvec
;   The problem I want to solve is how to calculate the
;   signal-to-noise ratio of the highest peak in the residual RM
;   spectrum if the covariance matrix of Q_ref and U_ref is not
;   diagonal and/or if the variances along the diagonal are not equal.
;   In general, confidence regions in the (Q,U) plane are nested ellipses 
;   that do not have to be aligned with the coordinate axes.
;   See appendix A in DS18.
;
;   The idea is that you can generate any distribution of data points
;   and their covariance matrix from a 2D distribution of points  
;   described by the covariance matrix [[1,0],[0,1]] by applying a  
;   rotation and a scaling. 
;   See Spruyt, 'A geometric interpretation of the covariance matrix', or
;   https://en.wikipedia.org/wiki/Whitening_transformation . 
;   What I'm doing here is known as PCA whitening.
;   If the scaling matrix is S, and the rotation matrix is R, then you
;   can transform each data point by calculating R * S times the data
;   point (using vector notation for the data point); the covariance 
;   matrix associated with these data then becomes R * S * S * R^-1.
;   As Spruyt explains, the matrix formed by calculating S * S is a 
;   diagonal matrix with the eigenvalues of the new covariance matrix 
;   along the diagonal. The normalized eigenvectors of the new
;   covariance matrix make up the columns in R.
;
;   Turning this around, I can calculate the eigenvalues and
;   eigenvector of the covariance matrix, and from them, the scaling
;   matrix S (which has the square roots of the eigenvalues along the 
;   diagonal) and the rotation matrix R (contains the normalized 
;   eigenvectors as columns). Then I use R to derotate the confidence  
;   ellipses (and the polarization vector LL of the highest peak in the 
;   residual RM  spectrum), to make the confidence ellipses line up with 
;   the coordinate axes. Finally, I rescale the Q and U axes using the 
;   matrix S so that the confidence ellipses become circles. The length  
;   of the scaled and derotated vector polvec is then the signal-to-noise 
;   ratio I want to calculate.

    common machine_precision, eps_machine

;   Calculate the covariance matrix from noise_q_arr and noise_u_arr
;   at a particular RM and \alpha = 0:
    cmatrix=Calc_noise_covariance(rm)
;
;   Since cmatrix is a 2x2 matrix, calculate its two eigenvalues and 
;   eigenvectors:
    Eigen22, cmatrix, val1=val1, vec1=vec1, val2=val2, vec2=vec2, eps_machine=eps_machine
;   Construct the scale and derotation matrices from the eigenvalues
;   and eigenvectors. First, normalize the eigenvectors:
    vec1n= vec1/sqrt(vec1[0]^2+vec1[1]^2)
    vec2n= vec2/sqrt(vec2[0]^2+vec2[1]^2)
;   The rotation matrix is equal to the matrix with normalized eigenvectors as its columns:
    rmatrix= [ [vec1n[0],vec2n[0]],[vec1n[1],vec2n[1]] ]  ; the rotation matrix
;   Test: does rmatrix times matrix produce a diagonal matrix, as it
;   should if rmatrix is indeed a rotation matrix? In that case
;   Inv(rmatrix) ## matrix rotates 'matrix' to the coordinate axes.
;    print, Inv22(rmatrix)##rmatrix
    smatrix= [ [sqrt(val1),0],[0,sqrt(val2)] ]  ; the scale matrix
;   Now calculate the matrix that you apply to the observed
;   polarization vector (see also Fig. 10 in Spruyt's notes):
    newvec= Inv22(rmatrix##smatrix)##polvec
    snr=sqrt(newvec[0]^2+newvec[1]^2)

    return, snr
END

FUNCTION FIND_MAX_CMP
;   Find the total number of model components you want to fit.
;   Note: max_cmp can be 0.

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    sel_tmp=where(model_type_arr ne 0, n_cmp)   
    return,n_cmp-1
END

FUNCTION CREATE_P_BIG, p, hat_vec_p_ref, max_cmp=max_cmp
;   Create an array 'p_big' which, for each source component, combines
;   information in 'p' with estimates of the Stokes Q and U
;   parameters at the reference frequency.

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    if n_elements(max_cmp) eq 0 then max_cmp=Find_max_cmp()

    p_big=dblarr((max_cmp+1)*n_param_max) 
    for ii=0,max_cmp do $
      p_big[ii*n_param_max : (ii+1)*n_param_max-1] = $
        [p[ii*(n_param_max-2) : (ii+1)*(n_param_max-2)-1], hat_vec_p_ref[2*ii], hat_vec_p_ref[2*ii+1]]

    return, p_big
END

FUNCTION CALC_LL, $
    stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_mod_tmp, stokes_u_mod_tmp, $
    n_param_tot=n_param_tot, chi2_red=chi2_red

;   Calculate the value of the log likelihood and the reduced chi
;   squared (if required).

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

;   Calculate the log likelihood:
    tmp= -0.5D * total(((stokes_q_arr_orig-stokes_q_mod_tmp)/noise_q_arr)^2,/double) $
         -0.5D * total(((stokes_u_arr_orig-stokes_u_mod_tmp)/noise_u_arr)^2,/double) 
;   Note that if you add 'eta' to 'tmp', you might need to change 'chi2_red' as well.
    logl= tmp $
          -total(alog(noise_q_arr),/double)-total(alog(noise_u_arr),/double) $
          -n_chan*alog(2*!dpi)

    if n_elements(chi2_red) eq 0 then return, logl
;   In this case, calculate also the reduced chi squared:
    chi2_red= -2 * tmp / (2l*n_chan-n_param_tot)  
;   Note that the number of observations is 2*n_chan  
    return, [logl, chi2_red]
END

FUNCTION SIMPLE_GAUSSIAN, p_arr, ampl_only=ampl_only
;   Return Stokes Q and U as a function of frequency for a point source or Gaussian RM distribution.
;
;   Note: p_arr has a length of n_param_max cells, and describes all 
;   parameters of a single source component.
;   p_arr can be a subarray of p_big.
;
;   p_arr[0] = RM, p_arr[1] = sigma_RM_intern (not used), p_arr[2] = \Delta_RM (not used),
;   p_arr[3] = sigma_RM_extern, p_arr[4] = \alpha, p_arr[5] = Q_ref, p_arr[6] = U_ref

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    if ~keyword_set(ampl_only) then begin
      arg= 2*p_arr[0]*lambda2_arr  ; 2*RM*lambda^2
      vec_P_mod = dcomplex(p_arr[5],p_arr[6])*(freq_arr/freq_ref)^p_arr[4]* $
        exp(-2*p_arr[3]^2*lambda2_arr^2)*dcomplex(cos(arg),sin(arg))
    endif $
    else begin
;     Consider only how the amplitude of the signal (normalised to
;     one) changes with frequency, as expressed by 'fnc_j' in section
;     2.2 of DS18.
      vec_P_mod=(freq_arr/freq_ref)^p_arr[4]*exp(-2*p_arr[3]^2*lambda2_arr^2)
    endelse

    return, vec_P_mod
END

FUNCTION SIMPLE_IFD, p_arr, ampl_only=ampl_only
;   Return Stokes Q and U as a function of frequency for an internal
;   Faraday dispersion model with an external Gaussian turbulent
;   screen, or simpler models. 
;   Also Burn slabs are calculated this way (use sigma_RM_intern=0).
;   'RM', the first parameter in the array p_arr, can be tricky to
;   interpret if the source emits at more than one RM. In that case
;   the parameter RM refers to the variables RM_0 or RM_c that are
;   used in section 2.1 in DS18.
;
;   Allow for sigma_RM_intern = 0, sigma_RM_extern = 0, and/or \Delta RM = 0.
;
;   Note: 
;   p_arr has a length of n_param_max cells, and describes all parameters of a single source component.
;   p_arr can be a subarray of p_big.
;
;   p_arr[0] = RM, p_arr[1]=sigma_RM_intern, p_arr[2]=\Delta_RM, p_arr[3]=sigma_RM_extern, 
;   p_arr[4] = alpha, p_arr[5] = Q_ref, p_arr[6] = U_ref

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    if ~keyword_set(ampl_only) then begin
      if ~(p_arr[1] eq 0D and p_arr[2] eq 0D) then begin
        arg= 2*p_arr[0]*lambda2_arr  ; 2*RM*lambda^2
;       ss = 2*sigma_RM_intern^2*lambda^4 - 2*i*\Delta_RM*lambda^2
        ss = dcomplex(2*p_arr[1]^2*lambda2_arr^2, -2*p_arr[2]*lambda2_arr)  
        vec_P_mod = dcomplex(p_arr[5],p_arr[6])*(freq_arr/freq_ref)^p_arr[4]* $
          Conj(ss)*(1-exp(-ss))/abs(ss)^2 * dcomplex(cos(arg),sin(arg)) * exp(-2*p_arr[3]^2*lambda2_arr^2)
      endif $ 
      else begin
;       If p_arr[1] eq 0D and p_arr[2] eq 0D then the variable 'ss' in the
;       denominator of \vec{P} becomes zero, causing numerical problems.
;       Avoid this by calling Simple_Gaussian in this case:
        vec_P_mod= Simple_Gaussian(p_arr)
      endelse
    endif $  ; if ~keyword_set(ampl_only)
    else begin
       if p_arr[1] eq 0D and p_arr[2] eq 0D then begin  
;       This means sigma_RM_intern eq 0 and \Delta_RM eq 0, hence call Simple_Gaussian:
        vec_P_mod= Simple_Gaussian(p_arr,/ampl_only)
      endif $
      else begin
        print,'Simple_IFD with /ampl_only does not work..'
        close,/all & stop
      endelse
    endelse
   
    return, vec_P_mod
END

FUNCTION FIND_STOKES_Q_U_MOD, fnc_name_arr, p_big, n_param_max, max_cmp=max_cmp
;   Model the Stokes Q and U frequency spectra, summing model
;   components between index 0 and 'max_cmp'.

    if n_elements(max_cmp) eq 0 then $
      Message,'Find_stokes_q_u_mod: please specify max_cmp'

    stokes_q_mod=0D & stokes_u_mod=0D  ; initialize
    for ii=0,max_cmp do begin
;     Predict the contribution by this source component to Stokes Q
;     and U as a function of frequency:
      vec_P_mod= Call_function(fnc_name_arr[ii], p_big[ii*n_param_max : (ii+1)*n_param_max-1])
;     Add this contribution vectorially to the previously calculated Stokes Q and U spectra:
      stokes_q_mod+= Real_part(vec_P_mod) & stokes_u_mod+= Imaginary(vec_P_mod)    
    endfor

    return, [stokes_q_mod,stokes_u_mod]
END

FUNCTION CALC_CF_SF_I, p_big, kk
;   Calculate the variables 'cf_i' and 'sf_i' from equation 11 in
;   DS18, which are used to determine the intrinsic polarization 
;   parameters Q_ref, U_ref of each source component.
;   'kk' is an index for the array model_type_arr.
;   Note that even though Calc_cf_sf_i is called using the array
;   p_big, which contains cells for holding Q_ref and U_ref, the
;   values inside these cells are not used by Calc_cf_sf_i,
;   therefore you can use 'p_big' safely in Calc_cf_sf_i.

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    index_kk = n_param_max*kk  ; index for the array p_big

    if (model_type_arr[kk] le 2) or $
       ((model_type_arr[kk] ge 3 and model_type_arr[kk] le 8) and $
        (p_big[index_kk+1] eq 0D and p_big[index_kk+2] eq 0D)) then begin  
;     p_big[index_kk+1] = sigma_RM_intern, p_big[index_kk+2] = delta_RM.
;     If these parameters are both zero then the denominator in the
;     IFD models equals zero. 
;     Avoid this by calling Simple_Gaussian:
      fnc= Simple_Gaussian(p_big[index_kk : index_kk + n_param_max-1], /ampl_only)
      arg= 2*p_big[index_kk]*lambda2_arr  ; rm=p_big[index_kk]
      cf_i= cos(arg)*fnc & sf_i= sin(arg)*fnc  ; Note: these are 1D arrays
    endif $
    else begin
;     In this case the model is for internal Faraday dispersion with
;     an external depolarizing Faraday screen, possibly with
;     sigma_RM_intern = 0, sigma_RM_extern = 0 and/or \Delta_RM = 0.
;     Note that you need to use 4 different functions (in two cases with a
;     non-zero \chi_0) to model the polarization vector from internal
;     Faraday dispersion.
;
;     Note: to improve numerical stability the case where both sigma_RM_intern=0
;     and \Delta_RM=0 (= turning IFD into an external DP screen or even a
;     point source) is passed to Simple_Gaussian to avoid the variable
;     'ss' in the denominator of \vec{P} becoming zero.
      rm= p_big[index_kk]
      sigma_RM_intern= p_big[index_kk+1]
      delta_RM= p_big[index_kk+2]
      sigma_RM_extern= p_big[index_kk+3]
      alpha= p_big[index_kk+4]
      denom= (2*sigma_RM_intern^2*lambda2_arr^2)^2 + (2*delta_RM*lambda2_arr)^2

      f1_i= 2*sigma_RM_intern^2*lambda2_arr^2/denom * (freq_arr/freq_ref)^alpha * exp(-2*sigma_RM_extern^2*lambda2_arr^2)
      arg= 2*rm*lambda2_arr
      cf1_i= cos(arg)*f1_i & sf1_i= sin(arg)*f1_i 

      f2_i= 2*delta_RM*lambda2_arr/denom * (freq_arr/freq_ref)^alpha * exp(-2*sigma_RM_extern^2*lambda2_arr^2)
      arg= 2*rm*lambda2_arr+!dpi/2
      cf2_i= cos(arg)*f2_i & sf2_i= sin(arg)*f2_i

      f3_i= -f1_i*exp(-2*sigma_RM_intern^2*lambda2_arr^2) 
      arg= 2*(rm+delta_RM)*lambda2_arr
      cf3_i= cos(arg)*f3_i & sf3_i= sin(arg)*f3_i

      f4_i= -f2_i*exp(-2*sigma_RM_intern^2*lambda2_arr^2) 
      arg= 2*(rm+delta_RM)*lambda2_arr+!dpi/2
      cf4_i= cos(arg)*f4_i & sf4_i= sin(arg)*f4_i

      cf_i = cf1_i + cf2_i + cf3_i + cf4_i
      sf_i = sf1_i + sf2_i + sf3_i + sf4_i
    endelse

    return, [[cf_i], [sf_i]]
END

FUNCTION FIND_MLE_Q_U_REF, $
    p_big, n_cmp=n_cmp, stokes_q_arr_orig=stokes_q_arr_orig, stokes_u_arr_orig=stokes_u_arr_orig

;   Find the maximum likelihood estimators for Q_ref and U_ref for a
;   given RM, sigma_RM (, ..), and spectral index alpha (specified in the
;   array 'p_big').

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common indices, cmp_index

    common flag_status, flag_status  ; stores error code

    common machine_precision, eps_machine

    if n_elements(n_cmp) eq 0 then n_cmp=cmp_index+1
;   here you need the number of components up to this cmp_index, not
;   the total number of all components in the model.
    data_vector=dblarr(2*n_cmp)  ; each source component specifies one Q_ref and U_ref

;   Calculate the elements of the matrix A_{ii,jj} on the left-hand side of
;   equation 11 in DS18 and the vector E_j from that paper.
;   Use indices ii and jj to cycle through all matrix elements.
;   Take care about the order of indices in equation 11 and how IDL
;   labels rows and columns. If there are two source components,
;   A_{1,2} is the element in the top-right of the matrix A in
;   equation 11 of DS18. Its IDL row index is 0 and column index
;   is 1, or, in IDL-speak, this is element "A[1,0]": exactly the
;   opposite order of indices compared to the notation in DS18.
;   (after compensating for the fact that IDL indices start at 0 and not 1).
    matrix=dblarr(2*n_cmp,2*n_cmp)  
    for ii=0,n_cmp-1 do begin
;     Calculate the variables 'cf1_i' and 'sf1_i': 
      tmp= Calc_cf_sf_i(p_big, ii)
;     Note: Q_ref and U_ref (set to zero in p_big) are not used in
;     Calc_cf_sf_i, therefore you can pass 'p_big' to that function.
      cf1_i=tmp[*,0] & sf1_i=tmp[*,1]

;     Because of symmetry you can start at jj=ii instead of jj=0 (see DS18).
      for jj=ii,n_cmp-1 do begin
;       Calculate the matrix aa = A[jj,ii] (IDL notation) = A_{ii+1,jj+1} 
;       (notation from the paper).
;       Use the matrix notation from DS18.

        if jj ne ii then begin
          tmp= Calc_cf_sf_i(p_big, jj)
          cf2_i=tmp[*,0] & sf2_i=tmp[*,1]
        endif $ 
        else begin
          cf2_i=cf1_i & sf2_i=sf1_i
        endelse

        bb= total( (cf1_i*cf2_i/noise_q_arr^2 + sf1_i*sf2_i/noise_u_arr^2) ,/double)

        cc1= total( (-cf1_i*sf2_i/noise_q_arr^2 + sf1_i*cf2_i/noise_u_arr^2) ,/double)

        cc2= total( (-cf2_i*sf1_i/noise_q_arr^2 + sf2_i*cf1_i/noise_u_arr^2) ,/double)        
;       Note that in general cc2 is not equal to -cc1.

        dd= total( (sf1_i*sf2_i/noise_q_arr^2 + cf1_i*cf2_i/noise_u_arr^2) ,/double)

        aa= [[bb, cc1],[cc2, dd]]  ; the matrix A[jj,ii] (IDL notation)

        matrix[2*jj : 2*jj+1, 2*ii : 2*ii+1] = aa
;       use symmetry to fill in cells below the diagonal of the matrix:
        if n_cmp gt 1 and ii ne jj then matrix[2*ii : 2*ii+1, 2*jj : 2*jj+1] = Transpose(aa)         
      endfor  ; for jj=ii,n_cmp-1
      vec_E = total(dcomplex(cf1_i,-sf1_i)*vec_p_freq_weighted,/double)
      data_vector[2*ii : 2*ii+1]= [Real_part(vec_E), Imaginary(vec_E)]  
;     'data_vector' is the vector on the right-hand side of the matrix equation 
;     to find (Q_ref,U_ref) of each source component (equation 8 in DS18).
    endfor  ; for ii=0,n_cmp-1

    n_pts=(2*long(n_cmp))^2  ; the number of elements in the matrix you want to invert
    sel=where(abs(matrix)/n_pts lt eps_machine, n_sel) 
    if n_sel gt 0 then matrix[sel]=0D
    if n_sel eq n_pts then return, [0D, 0D]  
;   Qref, Uref=0 means that the model component isn't there:
;   this avoids setting flag_status=-1000
;
;   Reduce the numbers in the matrix before inverting:
;   Use matrix^-1 = (matrix/N_pts)^-1/N_pts
    inv_matrix=Invert(matrix/n_pts,status,/double)/n_pts
    if status ne 0 then begin
      print,' Warning: matrix inversion unreliable..' 
      printf,99,'  '+Strjoin(roundoff(model_type_arr),' ')+' warning: matrix inversion unreliable..' 
      printf,99,matrix/n_pts
      flag_status=-1000
      sol_arr=dblarr(2,n_cmp)+flag_status
      if n_cmp gt 1 then return, sol_arr else return,sol_arr[*,0]
    endif
    sel=where(abs(inv_matrix) lt eps_machine, n_sel) 
    if n_sel gt 0 then inv_matrix[sel]=0D
    if flag_status eq -1000 then flag_status=0  
;   Sometimes flag_status eq -1000 with the first iteration of MPFIT:
;   if you don't reset flag_status here, flag_status will remain -1000
;   even if the matrix inversion was reliable.
    goto,skip_tests  ; comment out this line for debugging purposes.
    print,'------------'
    print,matrix 
    print,'111111111111'
    print,matrix/n_pts  ; this is the actual matrix that is being inverted
    print,'222222222222'
    print,inv_matrix*n_pts  ; this is the actual result from the matrix inversion
    print,'333333333333'
    tmp_matrix=matrix##inv_matrix
    sel_tmp=where(tmp_matrix lt eps_machine,n_sel_tmp)
    if n_sel_tmp gt 0 then tmp_matrix[sel_tmp]=0D
    print,tmp_matrix  ; should be equal to the identity matrix
    print,'============'
    print,''
    skip_tests:
    solution=inv_matrix##data_vector  
    solution=reform(solution[0,*]) 

    for ii=0, n_cmp-1 do $
       p_big[(ii+1)*n_param_max-2 : (ii+1)*n_param_max-1]=[solution[ii*2],solution[ii*2+1]] 
END

FUNCTION FIND_ERR_Q_U_REF, $
    stokes_q_arr_orig, stokes_u_arr_orig, p_big_cp, cmp_index=cmp_index, max_cmp=max_cmp, find_err_u=find_err_u
   
;   Calculate the error in Q_ref and U_ref of each model component.
;   'cmp_index' specifies for which source component you want to find the error in hat_q or hat_u.
;   'max_cmp' provides the index of the last source component you want to include in the model. 

;   Since p_big is modified in this function, work with a copy called 'p_big_cp'.

;   First, find a range in Stokes Q or U in such a way that you can fit a
;   parabola around the maximum in log likelihood with 9 data points.
;   The measurement errors in Q_ref and U_ref are based on this fit.

;   Note: check that stokes_q_mod and stokes_u_mod are calculated just
;   before Find_err_q_u_ref is called; otherwise the program might not 
;   calculate the correct log likelihood.

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod
;
    common param_noise_estimate, delta_logl_onesigma
;
    common flag_status, flag_status  ; stores error code
;
    logl_opt= Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_mod, stokes_u_mod)
    logl_onesigma=logl_opt-delta_logl_onesigma  ; contour level in log likelihood space

;   Find the error in the ML estimator for Q_ref or U_ref;
;   start with finding the low-end of the interval.
    stokes_index=cmp_index*n_param_max + 5 + keyword_set(find_err_u) 
;   = the index of Q_ref or U_ref in p_big of the component you're investigating
    stokes_opt= p_big_cp[stokes_index]  ; selects the value of Q_ref or U_ref
;   For the calculation of crude_err, see the covariance matrix in
;   equation 17 from Schnitzeler & Lee (2017) and the erratum for
;   that paper.
    alpha_opt= p_big_cp[stokes_index-1-keyword_set(find_err_u)]  ; selects the value of the spectral index \alpha.      
    if ~keyword_set(find_err_u) then $
      crude_err= Calc_crude_q_u_errors(noise_q_arr,alpha_opt) $
    else $
      crude_err= Calc_crude_q_u_errors(noise_u_arr,alpha_opt)
;   (see also Numerical Recipes in C, equation 15.6.4)
;   Assume that crude_err is a good approximation to the real error,
;   which means that the true error is not larger than 4*crude_err.
    stokes_range_min = stokes_opt - 4*crude_err
;   Probably crude_err is much closer to the true_error, with the
;   'while' loop below you'll quickly find a better value for stokes_range_min.       
;   Calculate the log likelihood of this q_range_min or u_range_min;
;   overwrite stokes_range_min at the position in p_big_cp used to store
;   Q_ref or U_ref of the model component you're fitting:
    p_big_cp[stokes_index]=stokes_range_min
    stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, p_big_cp, n_param_max, max_cmp=max_cmp)
    ll_stokes_range_min= $
      Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_u_mod[0 : n_chan-1], stokes_q_u_mod[n_chan : 2*n_chan-1])
;
    if ll_stokes_range_min gt logl_onesigma then $
;     q_range_min or u_range_min already falls inside the fit range of the
;     parabola, so there is no need to find a better q_range_min or u_range_min.
      goto,skip_refining_stokes_range_min

;   Find a better value for stokes_range_min:
    count=0  ; counter, keeps track of the number of iterations needed to find stokes_range_min.
    stokes_range_new= stokes_opt-(stokes_opt-stokes_range_min)/2D & count++
    ctu=1  ; boolean
    while ctu do begin
;     Refine q_range_min or u_range_min only while this parameter lies outside the
;     fit range of the parabola.

      p_big_cp[stokes_index]=stokes_range_new
      stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, p_big_cp, n_param_max, max_cmp=max_cmp)
      ll_stokes_range_new= $
        Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_u_mod[0 : n_chan-1], stokes_q_u_mod[n_chan : 2*n_chan-1])
      if count gt 20 then begin
        flag_status=-100
        return,[-flag_status,flag_status]  ; [err_q_min, err_q_max]
      endif
      if ll_stokes_range_new lt logl_onesigma then begin
        stokes_range_new= stokes_opt-(stokes_opt-stokes_range_new)/2D & count++ 
      endif $
      else begin
;       You've found a better value for stokes_range_min, so exit.
        stokes_range_min=stokes_range_new & ll_stokes_range_min=ll_stokes_range_new
        ctu=0 
      endelse
    endwhile  ; while ctu
    skip_refining_stokes_range_min:

;   Repeat for stokes_range_max:
    stokes_range_max = stokes_opt + 4*crude_err
    p_big_cp[stokes_index]=stokes_range_max
    stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, p_big_cp, n_param_max, max_cmp=max_cmp)
    ll_stokes_range_max= $
      Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_u_mod[0 : n_chan-1], stokes_q_u_mod[n_chan : 2*n_chan-1])

    if ll_stokes_range_max gt logl_onesigma then goto,skip_refining_stokes_range_max

    count=0  ; reset
    stokes_range_new= stokes_opt+(stokes_range_max-stokes_opt)/2D & count++
    ctu=1 
    while ctu do begin
      p_big_cp[stokes_index]=stokes_range_new
      stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, p_big_cp, n_param_max, max_cmp=max_cmp)
      ll_stokes_range_new= $
        Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_u_mod[0 : n_chan-1], stokes_q_u_mod[n_chan : 2*n_chan-1])
      if count gt 20 then begin
        flag_status=-100
        return,[-flag_status,flag_status]  ; [err_q_min, err_q_max]
      endif
      if ll_stokes_range_new lt logl_onesigma then begin
        stokes_range_new= stokes_opt+(stokes_range_new-stokes_opt)/2D & count++ 
      endif $
      else begin
        stokes_range_max=stokes_range_new & ll_stokes_range_max=ll_stokes_range_new
        ctu=0 
      endelse
    endwhile  ; while ctu
    skip_refining_stokes_range_max:

 
;   Find the error in Q_ref or U_ref by fitting a 1D parabola to the log
;   likelihood as a function of Q_ref or U_ref:
    stokes_arr=dblarr(9) & stokes_arr[0]=stokes_range_min  ; '9' = 9 sampling points
    ll_arr=dblarr(9) & ll_arr[0]=ll_stokes_range_min
    stokes_step=(stokes_range_max-stokes_range_min)/8D
    for ii=1,8 do begin
;     Note: ii=8 -> ss=stokes_range_max
      ss=stokes_range_min+ii*stokes_step
      stokes_arr[ii]=ss
      p_big_cp[stokes_index]=ss
      stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, p_big_cp, n_param_max, max_cmp=max_cmp)
      ll_arr[ii]= $
        Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_u_mod[0 : n_chan-1], stokes_q_u_mod[n_chan : 2*n_chan-1])
    endfor
    sel_tmp=where(ll_arr eq !VALUES.F_NAN, n_sel_tmp)
    if n_sel_tmp gt 0 then flag_status=-1
    sel_ss=where(ll_arr ge logl_opt-1D, n_sel_ss)
    if n_sel_ss eq 0 then flag_status=-2
    result=Fit_least_sq_solution(stokes_arr[sel_ss], ll_arr[sel_ss], /silent, /return_all_coefficients, /nointerrupt)
;   Note: you're allowed to call this function because the errors in
;   the y coordinates of the parabola fit are all zero.
    if finite(result[0],/nan) then flag_status=-3
    if flag_status ne 0 then return,[-flag_status,flag_status]

    stokes_low=result[1]-sqrt(-delta_logl_onesigma/result[2]) 
    stokes_high=result[1]+sqrt(-delta_logl_onesigma/result[2]) 
;   Convert from Q_min, Q_max (or U_min, U_max) to the measurement
;   uncertainty, which is the distance to Q_ref (or U_ref):
    err_stokes_low=-(result[1]-stokes_low) & err_stokes_high=stokes_high-result[1]
    err_arr=[err_stokes_low,err_stokes_high]

    sel_tmp=where(Finite(err_arr,/nan),n_sel_tmp)
    if n_sel_tmp eq 0 then return,err_arr $
    else begin
      flag_status=-4
      return,[-flag_status,flag_status]
    endelse
END

FUNCTION FIND_ERR_Q_U_REF_WRAPPER, stokes_q_arr_orig, stokes_u_arr_orig, p_big, cmp_index=cmp_index
;   Estimate the error in Q_ref and U_ref of each model component.
;   Note: run Find_err_q_u_ref only after you've finished
;   calculating the ML estimators for the last component in a source.
;   Don't pass cmp_index to this function using a common block
;   to avoid updating this parameter through call-by-reference (see
;   the double for loop below).

    n_cmp=cmp_index+1  ; include model components up to and including 'cmp_index'

    err_arr=dblarr(4,n_cmp)  ; two errors per Stokes per model component
;   Note: Find_err_q_u_ref finds the error in either Q_ref or U_ref
;   for one source component at a time. Therefore you need the
;   following double for loop to calculate the measurement
;   uncertainties in these parameters for all source components:
    for ii=0,cmp_index do $
      for jj=0,1 do begin
;       p_big will be modified by Find_err_q_u_ref, therefore work
;       with a copy of this array:
        p_big_cp=p_big
        err_arr[2*jj : 2*jj+1,ii]= $
          Find_err_q_u_ref(stokes_q_arr_orig, stokes_u_arr_orig, p_big_cp, $
                           cmp_index=ii, find_err_u=jj, max_cmp=cmp_index)

      endfor  ; for jj=0,1
    if n_cmp eq 1 then err_arr=err_arr[*,0]

    return, err_arr
END

FUNCTION MODEL_WRAPPER, x, p
    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common indices, cmp_index
;   keep this common block for passing 'cmp_index' to Model_wrapper. 
;   MPFITFUN calls Model_wrapper, passing only 'x' and 'p' to this function.

    common machine_precision, eps_machine

    common p_big, p_big

;   Note: we determine the ML estimators for Q_ref and U_ref using
;   linear algebra, using the non-linear parameters as input.
;   Therefore, MPFITFUN provides the function Model_wrapper with
;   guesses for the non-linear model parameters in the array 'p', and
;   Model_wrapper calculates the ML estimators for Q_ref and U_ref
;   based on 'p'.
;
;   'x' is a dummy variable,
;   'p' is an array which contains only(!) the non-linear parameters of each model component. 

;   Find the number of model components you want to fit:
    n_cmp=cmp_index+1
;   here you need the number of components up to this cmp_index, not
;   the total number of all components in the model.

;   If you're fitting an IFD-based model, check if calculating sigma_RM_intern^2 + delta_RM^4 leads to
;   problems with the numerical accuracy:
    sel=where(model_type_arr ge 3, n_sel)
    if n_sel gt 0 then begin
      index=0
      for ii=0,cmp_index do begin
        if model_type_arr[ii] ge 3 then begin
          min_l2=min(lambda2_arr)
          ss_min = dcomplex(2*p[index+1]^2*min_l2^2, -2*p[index+2]*min_l2) 
;         This is the factor 'ss' in the denominator of an IFD model
;         (see the function Simple_IFD above), evaluated at the
;         smallest value of lambda squared of the observations. 
          if abs(ss_min)^2 lt 10*eps_machine then p[index+1 : index+2] = 0D
        endif
        index+=n_param_max-2
      endfor
    endif

;   Note: MPFITFUN only optimizes the non-linear model parameters,
;   therefore the array 'p' does not contain cells to store Q_ref and
;   U_ref. However, Find_mle_q_u_ref does expect an array with these
;   cells, even though the values for Q_ref and U_ref are currently
;   unknown. To solve this paradox, create an empty array and combine
;   this with 'p' into the array 'p_big', which is passed to
;   Find_mle_q_u_ref:
    p_big=Create_p_big(p, dblarr(2,n_cmp), max_cmp=cmp_index)  
    tmp=Find_mle_q_u_ref(p_big, n_cmp=n_cmp)
;   Find_mle_q_u_ref fills the cells in p_big which are supposed to
;   contain Q_ref and U_ref.
;
;   Model the Stokes Q and U frequency spectra, summing all model
;   components up to and including 'max_cmp':
    stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, p_big, n_param_max, max_cmp=cmp_index)
    stokes_q_mod=stokes_q_u_mod[0 : n_chan-1] & stokes_u_mod=stokes_q_u_mod[n_chan : 2*n_chan-1]

    return, [stokes_q_mod,stokes_u_mod]
END

PRO FIT_MODEL, stokes_q_arr_orig, stokes_u_arr_orig
;   Fit a model to the Stokes Q and U data for all model components
;   from zero up to and including cmp_index.

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common rm_alpha_grid, rm_min, rm_max, rm_step, alpha_min, alpha_max, alpha_step

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common indices, cmp_index

    common flag_status, flag_status  ; stores error code

    common p_big, p_big

;   Run RM synthesis on the residual Stokes Q,U frequency spectrum to
;   find a crude RM that can be used to fit the next model component:
    vec_p_curlyR = Calc_rms_freq(stokes_q_arr, stokes_u_arr, freq_arr, freq_step, $
                                 rm_min, rm_max, rm_step, /use_old_formalism)
;   'Calc_rms_freq' changes the dimensions of the data vectors; undo this:
    stokes_q_arr=stokes_q_arr[0:n_chan-1] & stokes_u_arr=stokes_u_arr[0:n_chan-1]
    crude_pi_max = max(abs(vec_p_curlyR), index) 
    crude_rm= index*rm_step + rm_min
;   Use this RM and alpha=-0.7d as starting points for the Levenberg-Marquardt
;   algorithm.
;   If a source is not a point source: set the initial guesses for
;   additional model components equal to 0.

;   Skip adding components if the highest peak in the residual RM spectrum
;   is only a few times the measurement uncertainty in Stokes Q and U:
;   (see also section 6.3 in DS18)
;    pi_cutoff= 1.5*sigma_q_rms  ; assuming sigma_q_rms = sigma_u_rms
;   Drop the requirement that sigma_q_rms = sigma_u_rms:
    snr_cutoff=1.5  ; only select peaks in the residual RM spectrum
                    ; with a signal-to-noise ratio >= pi_cutoff.
    tmp1=Real_part(vec_p_curlyR[index]) & tmp2=Imaginary(vec_p_curlyR[index])
    snr=Calc_snr(crude_rm, [tmp1,tmp2])
;
;   Note: because of noise bias, the polarized flux density of the
;   highest peak in the RM spectrum follows a more complex
;   distribution than a Rayleigh distribution: see Hales et
;   al. (2012), figure 2 in DS & Lee (2017). Therefore, one can choose
;   a larger value for pi_cutoff to avoid cleaning peaks in the RM
;   spectrum that are due to noise. By choosing a small value for
;   pi_cutoff we include as many model components as possible, and
;   rely on the AIC and BIC to select the model that best describes
;   the data. This allows us to investigate if the AIC and BIC always
;   select the injected model, and how they penalize models with
;   different numbers of parameters.
    if snr lt snr_cutoff then begin
      stokes_q_mod=dblarr(n_chan) & stokes_u_mod=dblarr(n_chan)
;     If stokes_q_mod and stokes_u_mod only contain zeroes (this way
;     Calc_metrics assigns this model a detection significance of
;     zero and a reduced chi squared of the fit >> 1), p_big, model_param,  
;     and model_err must contain zeroes for consistency:
      n_param_tmp= (cmp_index+1)*n_param_max  ; the total number of parameters of all components.
      p_big=dblarr(n_param_tmp)  
      model_param*=0 & model_err*=0
;
      flag_status = 111  
;     The name for this flag status refers to a list of point source
;     components, if this flag is set then one component has been
;     added too many.
      return
    endif

;   Define the starting point for the Levenberg-Marquardt algorithm, and which parameters should be kept fixed.
    model_name='Model_wrapper'  
;   Describe the starting point for the Levenberg-Marquardt algorithm.
;   First, calculate the total number of parameters (not counting
;   Q_ref and U_ref) for all source components combined:
    n_param_incl_this_cmp= (cmp_index+1)*(n_param_max-2)
;   '-2': Q_ref and U_ref are not specified in start_point.
    start_point=dblarr(n_param_incl_this_cmp)
    start_point[n_param_incl_this_cmp - (n_param_max-2) : n_param_incl_this_cmp-1]= $
      (fwhm_rmsf gt 0) ? $
      [crude_rm, replicate((0.83D*fwhm_rmsf/4), (n_param_max-2)-2), -0.7D] : $
      [crude_rm, replicate((0.83D*abs(fwhm_rmsf)/100), (n_param_max-2)-2), -0.7D]  
;   The starting guess for sigma_RM, delta_RM, ... is not 0
;   but a fraction of the resolution of the RMSF.
;   (this way the algorithm works better)
;   '(n_param_max-2)' since you only specify starting values for the
;   non-linear model parameters, an extra '-2' because you specified
;   crude_rm and alpha=-0.7 separately.
;   Note variables that are not used in the fitting are set to zero
;   in the for loop a few lines lower.
;   If Lookup_test_param sets fwhm_rmsf to -1 you force the initial
;   parameters for sigma_RM_internal, delta_RM, and sigma_RM_external
;   to zero.
    if cmp_index gt 0 then begin
;     Store also the components fitted in the parent model in start_point
      for ii=0,cmp_index-1 do $
        start_point[ii*(n_param_max-2) : (ii+1)*(n_param_max-2)-1]= model_param[0: (n_param_max-2)-1, ii]
;       '-2': you don't specify Q_ref and U_ref in start_point.
    endif

;   Fix parameters that are not used in the fitting, limit certain
;   other parameters:
    par_properties= replicate({fixed:0,limited:[0,0],limits:[0D,0D]}, n_param_incl_this_cmp)
    for ii=0,cmp_index do begin
      Lookup_model,model_type_arr[ii],param_used=param_used
      for jj=1,3 do begin
;       Note: jj=0 (RM foreground), 5,6 (Q_ref,U_ref) are not
;       constrained. jj=4 (spectral index alpha) is constrained 
;       after this for loop.
        if param_used[jj] eq 1 then begin
;         Limit the fit range for sigma_RM_intern, \Delta RM, and
;         sigma_RM_extern to values >= 0
          par_properties[ii*(n_param_max-2)+jj].limited[0]=1
;         '-2': Q_ref and U_ref are not specified in start_point
          par_properties[ii*(n_param_max-2)+jj].limits=0D
        endif $  ; if param_used[jj] eq 1
        else begin
;         Parameter is not used in the model, therefore fix it to zero
;         in the fitting:
          par_properties[ii*(n_param_max-2)+jj].fixed=1
          start_point[ii*(n_param_max-2)+jj]=0D
        endelse
      endfor  ; for jj=1,3
;     Only allow spectral index values between -6 and 3:
      par_properties[ii*(n_param_max-2)+4].limited[0:1]=1
      par_properties[ii*(n_param_max-2)+4].limits=[-6D,3D]
    endfor  ; for ii=0,cmp_index

    tmp_result=Mpfitfun(model_name,indgen(2*n_chan),[stokes_q_arr_orig,stokes_u_arr_orig], $
                               [noise_q_arr,noise_u_arr],start_point, parinfo=par_properties, $
                               covar=covar, perror=perror, maxiter=1000, $
                               status=status, errmsg=errmsg, silent=1)

    fit_result=p_big  ; tmp_result contains only the non-linear
                      ; parameters, fit_result includes all parameters.
                      ; Note: 'p_big' is passed back through the common block 'p_big'
    if status le 0 then Message,errmsg
;   If status le 0 then the program stops, so you don't need to assign
;   a special value to flag_status in that case.
    if status ne 1 then begin
      print,'Status: '+roundoff(status)
      printf,99,'  '+Strjoin(roundoff(model_type_arr),' ')+' status: '+roundoff(status)
    endif
    if status eq 5 then flag_status=-2000  ; this exit status from MPFIT can indicate that MPFIT failed to converge.
;    print,'Result: ',fit_result
;    print,' Covariance matrix: ' & print,covar

;   Find the model prediction for Stokes Q and U as a function of frequency:
    stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr, fit_result, n_param_max, max_cmp=cmp_index)
    stokes_q_mod=stokes_q_u_mod[0 : n_chan-1] & stokes_u_mod=stokes_q_u_mod[n_chan : 2*n_chan-1]

;   Find the errors for Q_ref and U_ref:
;   (Q_ref and U_ref themselves are returned in the array p_big by Fit_model)
    err_q_u_ref= $
      Find_err_q_u_ref_wrapper(stokes_q_arr_orig, stokes_u_arr_orig, fit_result, cmp_index=cmp_index)
;   Note: since this function depends on stokes_q_mod and
;   stokes_u_mod, check that you use up-to-date versions of these
;   arrays.
    for ii=0,cmp_index do begin
;     Note: the correct ML estimators for Q_ref and U_ref are stored
;     in 'fit_result', but the uncertainties of these parameters are not
;     in 'perror' but in 'err_q_u_ref'.
      model_param[*,ii]=fit_result[ii*n_param_max : (ii+1)*n_param_max-1]
      model_err[*,ii]= [perror[ii*(n_param_max-2) : (ii+1)*(n_param_max-2)-1], err_q_u_ref[0:3,ii]]
    endfor
END

FUNCTION CALC_METRICS, n_param_tot, stokes_q_arr_orig, stokes_u_arr_orig

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common flag_status, flag_status  ; stores error code

;   Calculate the log likelihood ratio, AIC (Akaike information criterion),  
;   BIC (Bayesian information criterion), the reduced chi squared, and
;   the maximum reduced chi squared:
    tmp= Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_mod, stokes_u_mod, n_param_tot=n_param_tot, /chi2_red)
    logl_opt= tmp[0] 
    chi2_red= tmp[1]
    logl0= Calc_ll(stokes_q_arr_orig, stokes_u_arr_orig, stokes_q_mod*0, stokes_u_mod*0)
    loglr= logl_opt-logl0  ; the log likelihood ratio
    delta_n_param=n_param_tot  ; the number of parameters used to
                               ; calculate log \Lambda - the number of
                               ; parameters used to calculate log \Lambda_0
    chi2_red_max=Find_chi2_red_max(n_param_tot=n_param_tot)

;   Note that the sample size (number of observations) used to calculate AIC and BIC is 2*n_chan.
    AIC= -2*logl_opt + 2*delta_n_param 
    BIC= -2*logl_opt + delta_n_param*alog(2*n_chan)

;   Is this a significant detection?
    prob=erf(5/sqrt(2D))  ; probability that |X| <= 5*sigma for a Gaussian distribution with mean 0 and sigma 1.
;   Calculate the value of the chi-squared distribution that produces the same significance:
;    level_chisq=Chisqr_cvf(1D -prob, delta_n_param) 
;   I used this method for the Monte Carlo simulations that I
;   described and analysed in DS18. However, 'chisqr_cvf'
;   hasn't been written yet for GDL. Use my own version:
    level_chisq=Chisqr_cvf_ds(prob, delta_n_param)
;   More accurate results can be obtained with Mathematica, if required.
    det_significance= (2*loglr gt level_chisq) ? 5 : -1

    return, [loglr, AIC, BIC, chi2_red, chi2_red_max, det_significance, flag_status]
END

PRO PLOTIT, $
    stokes_q_arr_orig, stokes_u_arr_orig, real_data=real_data, src=src, stokes_i_arr=stokes_i_arr, $
    noise_i_arr=noise_i_arr, chi2_red=chi2_red

;   Plot the data as a function of frequency, together with the
;   properties of the best-fitting model, as selected by the Bayesian
;   Information Criterion.
;   If the keyword 'real_data' is used, then Plotit will draw four
;   panels, from top to bottom:
;   (1) Measured Stokes I flux densities as a function of frequency,
;       together with 1-sigma uncertainties in these measurements.
;       Also shown are the best-fitting power-law model (found by
;       maximum likelihood). 
;       Two numbers are displayed in the top-right corner of this
;       panel: the first is the spectral index of this model,
;       \alpha_I, the second is the reduced chi squared of the fit of
;       the power-law to the Stokes I measurements.
;   (2) Stokes Q and U as a function of frequency (solid line and
;       dashed line, respectively). Overplotted is the best-fitting
;       model, as selected by the Bayesian Information Criterion. The
;       model prediction for Stokes Q is shown in blue, for Stokes U
;       it is shown in red.
;   (3) The residuals of Stokes Q and U as a function of frequency,
;       which shows the difference between the measurements and the
;       best-fitting model that was shown in the second panel.
;       The residuals in Stokes Q are shown in blue, the residuals in
;       Stokes U in red (the same colour coding as using in the second panel).
;       Error bars indicate 1-sigma errors.
;   (4) A summary of the model that fits the data best.
;       This panel shows how much polarized flux density is emitted at
;       each RM, not only for point sources, but also for rectangular
;       or Gaussian distributions of polarized flux density as a
;       function of RM (more complex models cannot be plotted yet). 
;       Note that the program does not add the polarization vectors of 
;       sources that overlap in RM (for example, two Burn slabs).
;       Also shown are two columns with numbers, together with 1-sigma
;       errors (where applicable).
;       In the column on the left, these numbers are:
;       * the Stokes I flux density at the reference frequency, which
;         was fitted to the Stokes I data in the panel at the top.
;       * the polarized flux density of the brightest component that
;         was fitted to the polarization data. The power law that we
;         fit to the data describes the intrinsic emission, before
;         this emission is depolarized. The number mentioned here
;         refers to the (intrinsic) polarized flux density of the
;         source at the chosen reference frequency.
;       * the polarization percentage of the brightest polarized
;         source component, again at the reference frequency. This is 
;         the ratio between the first two flux densities that are listed
;         in this column.
;       * the signal-to-noise ratio of the brightest polarized source
;         component, again at the reference frequency. This number is
;         calculated following the method outlined in appendix A of
;         DS18.
;       In the column on the right, either three or four numbers are
;       shown. These are:
;       * the RM of the brightest polarized source component. If not a
;         point source in RM but a different model was fitted, then the
;         number shown here is the RM_0 or RM_c used in sections 2.1 and
;         2.2 from DS18.
;       * the spectral index that describes the power law polarized
;         flux density spectrum of the brightest polarized source
;         component, \alpha_{L,1} (this emission is assumed to follow
;         a power law before it gets depolarized).
;       * the reduced chi squared of the model fitted to Stokes Q and U.
;       * the ratio between the polarized flux densities of the
;         brightest and second brightest component, again at the
;         reference frequency. This ratio gives you an idea if a
;         single component can account for the observed spectra of
;         Stokes Q and U, or if other components are important as well.
;         (this number is not shown if the best-fitting model
;          consists only of a single component)
;   If the keyword 'real_data' is not set, then Plotit shows the
;   data as a function of frequency. Stokes Q is shown with a solid
;   line, Stokes U with a dashed line. The coloured lines are the
;   predictions for Stokes Q and U by the best-fitting model: blue
;   shows the model prediction for Stokes Q, red for Stokes U.
;   If you are running a Monte Carlo simulation, then two additional
;   dashed lines are plotted. These are the predictions by the
;   injected (simulated) model for Stokes Q and U, before noise was
;   added. Blue is again for Stokes Q, red for Stokes U.
;
;   Note that the components listed with subscripts '1' and '2' in the 
;   PostScript file are the brightest and second brightest polarized
;   sources that were fitted to the data. Plotit sorts the list of
;   model components internally to find out which two polarized
;   sources are brightest; in general, in fs.pro the order of model 
;   components reflects the order in which they were fitted to the data.

    common input_model, model_type_arr_input, fnc_name_arr_input, p_big_input1, p_big_input2

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common report, report_i_param_err, report_model_param, report_model_err, report_snr, $
                   report_chi2_red, report_ind

    charsize= (!d.name ne 'PS') ? 2 : 1.2
    cs_mult= 0.8  ; the CharSize multiplier
    lt_mult= (!d.name ne 'PS') ? 1 : 3  ; the Line Thickness multiplier
    symsize= (!d.name ne 'PS') ? 1 : 0.6
    style= (!d.name ne 'PS') ? 1 : 2
    freq_range= Find_plot_range(freq_arr)

    loadct,39,/silent

    if keyword_set(real_data) then begin
      !p.multi=[0,1,4]

      yr=Find_plot_range([[stokes_i_arr-noise_i_arr],[stokes_i_arr+noise_i_arr]])
      plot,freq_arr,stokes_i_arr,psym=8,xr=freq_range,/xs,yr=yr,/ys,/Nodata, $
        xtitle='Frequency (MHz)',ytitle='Stokes I (mJy)', $
        title=((n_elements(src) gt 0) ? src : ''),charsize=charsize
      oplot,!x.crange,[0,0],linestyle=style,thick=1*lt_mult

      vsym,20,fill=(!d.name ne 'PS')
;     I learned the following speed-up trick from the GDL
;     implementation of ploterr.pro (by Gilles Duvert and Alain Coulais)
      null=replicate(!values.d_nan,n_chan)
      x_new=reform(transpose([[freq_arr],[freq_arr],[null]]),3*n_chan)
      oplot, freq_arr, stokes_i_arr, psym=8, symsize=symsize
      y_new=reform(transpose([[stokes_i_arr-noise_i_arr],[stokes_i_arr+noise_i_arr],[null]]),3*n_chan)
      plots,x_new,y_new,thick=1*lt_mult

;     Fit a power law to the Stokes I data:
      fit_result=Fit_powerlaw(freq_arr,stokes_i_arr,noise_i_arr,x_ref=freq_ref)
      stokes_i_ref=fit_result[0] & stokes_i_ref_err=fit_result[1]
      oplot,freq_arr,fit_result[0]*(freq_arr/freq_ref)^fit_result[2],thick=2*lt_mult,color=245
      report_i_param_err=dblarr(8)  
;     This array stores:
;     (0) flux density at the reference frequency 
;     (1) the associated error
;     (2) the spectral index 
;     (3) the associated error
;     (4) the reference frequency 
;     (5) the reduced chi squared of the fit
;     (6,7) the polarization percentages of the brightest and second
;           brightest polarized source at the reference frequency (filled in later)
      report_i_param_err[0:3]=fit_result[0:3]
      report_i_param_err[4]=freq_ref
      report_i_param_err[5]=fit_result[4]

      xr=!x.crange & yr=!y.crange
      xpos= 0.97*(xr[1]-xr[0])+xr[0]  ; check also the positions that are used in the fourth panel
      ypos= [0.85,0.72]*(yr[1]-yr[0])+yr[0]
      str='!4a!3!dI!n : '+roundoff(fit_result[2],2)+' !9+!3 '+roundoff(fit_result[3],2)
      xyouts,xpos,ypos[0],str,charsize=cs_mult*charsize,align=1
      str='!4v!3!u2!n!dred!n : '+roundoff(fit_result[4],2)
      xyouts,xpos,ypos[1],str,charsize=cs_mult*charsize,align=1
   endif

    yr= Find_plot_range([stokes_q_arr_orig,stokes_u_arr_orig,stokes_q_mod,stokes_u_mod])
    plot, freq_arr, stokes_q_arr_orig, xr=freq_range, /xs, yr=yr, /ys, /Nodata, $
      xtitle='Frequency (MHz)',ytitle='Stokes Q,U (mJy)', charsize=charsize, symsize=symsize
    oplot,!x.crange,[0,0],linestyle=style,thick=1*lt_mult
    oplot, freq_arr, stokes_q_arr_orig
    oplot, freq_arr, stokes_u_arr_orig, linestyle=2

;   Overplot the model fit:
    oplot,freq_arr,stokes_q_mod,color=80,thick=3*lt_mult  ; Blue = Stokes Q
    oplot,freq_arr,stokes_u_mod,color=245,thick=3*lt_mult
    if !d.name eq 'PS' then begin
      xpos_tmp=0.8*(!x.crange[1]-!x.crange[0])+!x.crange[0]
      ypos_tmp=0.95*(!y.crange[1]-!y.crange[0])+!y.crange[0]
      legend_pos=[xpos_tmp,ypos_tmp]
      Legend,['Q','U'],charthick=2,psym=[8,8],color=[80,245],/top_legend,/horizontal,pos=legend_pos,box=0
    endif

    if keyword_set(real_data) then begin
;     Plot the residuals in Stokes Q and U as a function of frequency:
      q_res= stokes_q_arr_orig-stokes_q_mod
      u_res= stokes_u_arr_orig-stokes_u_mod
      yr=Find_plot_range([q_res-noise_q_arr,q_res+noise_q_arr,u_res-noise_u_arr,u_res+noise_u_arr])
      plot, freq_arr, q_res, xr=freq_range, /xs, yr=yr,/ys, /Nodata, $
        xtitle='Frequency (MHz)',ytitle='Stokes Q,U residuals (mJy)', charsize=charsize
      oplot,!x.crange,[0,0],linestyle=style,thick=1*lt_mult
      if !version.os eq 'linux' then begin
;       You're running IDL
        vsym, 20, /fill, color=80
        oplot, freq_arr, q_res, psym=8, symsize=symsize
      endif $
      else begin
;       You're running GDL
        vsym,20,/fill
        oplot, freq_arr, q_res, psym=8, color=80, symsize=symsize
      endelse
      y_new=reform(transpose([[q_res-noise_q_arr],[q_res+noise_q_arr],[null]]),3*n_chan)
      plots,x_new,y_new,color=80,thick=1*lt_mult
      if !version.os eq 'linux' then begin
;       You're running IDL
        vsym, 20, /fill, color=245
        oplot, freq_arr, u_res, psym=8, symsize=symsize
      endif $
      else $
;       You're running GDL
        oplot, freq_arr, u_res, psym=8, color=245, symsize=symsize
      y_new=reform(transpose([[u_res-noise_u_arr],[u_res+noise_u_arr],[null]]),3*n_chan)
      plots,x_new,y_new,color=245,thick=1*lt_mult

;     Plot the polarized flux density of each source component as a function of RM:
      sel=where(model_err[0,*] ne 0,n_sel)  
;     this 'sel' removes rows from model_err that were not fitted,
;     selecting only components that were fitted.
      xr=[floor(min(model_param[0,sel],max=tmp_max)),tmp_max]+[-1,1]*5
      xr=Find_plot_range(xr)
      pi_arr=sqrt(model_param[5,*]^2+model_param[6,*]^2)
;     Sort pi_arr in order of descreasing polarized flux density:
      ind_sorted=reverse(bsort(pi_arr))
      yr=[0,1.2*pi_arr[ind_sorted[0]]]

      plot, xr, yr, xr=xr, /xs, yr=yr, /ys, /Nodata, $
        xtitle='RM (rad m!u-2!n)',ytitle='Flux density (mJy)', charsize=charsize
    for ii=0,n_sel-1 do $
;     At the moment only point sources, Gaussian, and rectangular distributions
;     of polarized flux density as a function of RM can be plotted.
      if model_param[1,sel[ii]] eq 0 then begin
        if model_param[1,sel[ii]] eq 0 and model_param[2,sel[ii]] eq 0 then $
;         Plot a point source in RM:
          oplot,[1,1]*model_param[0,sel[ii]],[0,pi_arr[sel[ii]]],thick=2*lt_mult,color=80
        if model_param[1,sel[ii]] ne 0 and model_param[2,sel[ii]] eq 0 then begin
;         Plot the Gaussian |L(RM)| distribution.
;         First: find out the scale
          delta_x=!x.crange[1]-!x.crange[0]
          x_step=delta_x/100.
          x_arr_tmp=findgen(4*model_param[3,sel[ii]]/x_step)*x_step
          y_arr_tmp=exp(-(x_arr_tmp/model_param[3,sel[ii]])^2/2.)
          oplot,model_param[0,sel[ii]]-x_arr_tmp,pi_arr[sel[ii]]*y_arr_tmp,thick=2*lt_mult,color=80
          oplot,model_param[0,sel[ii]]+x_arr_tmp,pi_arr[sel[ii]]*y_arr_tmp,thick=2*lt_mult,color=80
        endif  
        if model_param[1,sel[ii]] eq 0 and model_param[2,sel[ii]] ne 0 then begin
;         For the Burn slab, convert RM_0 and Delta_RM back to RM_min and RM_max:          
          rm_min= model_param[0,sel[ii]]-model_param[2,sel[ii]]/2.
          rm_max= rm_min + model_param[2,sel[ii]]
          oplot,[rm_min,rm_min,rm_max,rm_max],[0,pi_arr[sel[ii]],pi_arr[sel[ii]],0],thick=2*lt_mult,color=80
        endif
      endif
;     In this for loop it doesn't matter if the entries in model_param
;     have been sorted in order of descending polarized flux density.
      axis, xaxis=0, xr=xr, /xs, charsize=charsize

      xpos= [0.03,0.97]*(xr[1]-xr[0])+xr[0]
      ypos= [0.85,0.72,0.59,0.46,0.33]*(yr[1]-yr[0])+yr[0]
      str='I!dref!n : '+roundoff(stokes_i_ref)+' !9+!3 '+roundoff(stokes_i_ref_err)+' mJy'
      xyouts,xpos[0],ypos[0],str,charsize=cs_mult*charsize,align=0
      str='L!d1,ref!n : '+roundoff(pi_arr[ind_sorted[0]],1)+' !9+!3 '+ $
        roundoff((-model_err[5,ind_sorted[0]]+model_err[6,ind_sorted[0]])/2d,1)+' mJy'
      xyouts,xpos[0],ypos[1],str,charsize=cs_mult*charsize,align=0
      tmp=100d * pi_arr[ind_sorted[0]]/stokes_i_ref
      str='L!d1,ref!n/I!dref!n : '+roundoff(tmp,1)+' %'
      xyouts,xpos[0],ypos[2],str,charsize=cs_mult*charsize,align=0
      report_i_param_err[6]=tmp
;     calculate the signal-to-noise ratio using the full covariance matrix:
      polvec_tmp=[model_param[5,ind_sorted[0]],model_param[6,ind_sorted[0]]]
      snr1= Calc_snr(model_param[0,ind_sorted[0]], polvec_tmp)
      str='L!d1,ref!n/err : '+roundoff(snr1,1)
      xyouts,xpos[0],ypos[3],str,charsize=cs_mult*charsize,align=0

      str='RM!d1!n : '+roundoff(model_param[0,ind_sorted[0]],1)+ $
          ' !9+!3 '+roundoff(model_err[0,ind_sorted[0]],1)+' rad!u-2!n'
      xyouts,xpos[1],ypos[0],str,charsize=cs_mult*charsize,align=1
      str='!4a!3!dL,1!n : '+roundoff(model_param[4,ind_sorted[0]],2)+ $
          ' !9+!3 '+roundoff(model_err[4,ind_sorted[0]],2)
      xyouts,xpos[1],ypos[1],str,charsize=cs_mult*charsize,align=1
      str='!4v!3!u2!n!dred!n : '+roundoff(chi2_red,2)
      xyouts,xpos[1],ypos[2],str,charsize=cs_mult*charsize,align=1
;     Note: if a component that is more complex than a point
;           source was fitted then you're only comparing the peak
;           amplitude of that component with the peak amplitude of the 
;           first model component. You're not comparing integrated flux densities.
      if n_sel ge 2 then begin
        str='L!d1,ref!n/L!d2,ref!n : '+roundoff(pi_arr[ind_sorted[0]]/pi_arr[ind_sorted[1]],1)
        xyouts,xpos[1],ypos[3],str,charsize=cs_mult*charsize,align=1
;       Calculate the polarization percentage of the second brightest
;       source component:
        report_i_param_err[7]=100d * pi_arr[ind_sorted[1]]/stokes_i_ref
;       Also calculate the signal-to-noise ratio for the second brightest peak:
        polvec_tmp=[model_param[5,ind_sorted[1]],model_param[6,ind_sorted[1]]]
        snr2= Calc_snr(model_param[0,ind_sorted[1]], polvec_tmp)
      endif $
      else begin
        report_i_param_err[7]=-1
        snr2=-1
      endelse
      report_snr=[snr1,snr2]

      !p.multi=0
    endif $
    else begin
;     Overplot the model you used to generate the mock data set:
      sel_tmp=where(model_type_arr_input ne 0, n_cmp_inj_max)
      p_big= (n_cmp_inj_max eq 2) ? [p_big_input1,p_big_input2] : p_big_input1
      stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr_input, p_big, n_param_max, max_cmp=n_cmp_inj_max-1)
;
      oplot,freq_arr,stokes_q_u_mod[0 : n_chan-1],color=80,thick=3,linestyle=2
      oplot,freq_arr,stokes_q_u_mod[n_chan : 2*n_chan-1],color=245,thick=3,linestyle=2
    endelse

    loadct,0,/silent
END

FUNCTION FIND_MODEL_STR
;   Combine information on this model into a single string.

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    max_cmp=Find_max_cmp()

    for ii=0,max_cmp do begin
;     Collect information for this model component:
      model_str=''
      for jj=0,n_param_max-3 do begin  ; loop over all parameters except Q_ref and U_ref
        model_str+= roundoff(model_param[jj,ii],5)+' '+roundoff(model_err[jj,ii],5)+' '
      endfor
      for jj=0,1 do begin
;       Special loop to store Q_ref and its two associated errors (same for U_ref).
        model_str+= roundoff(model_param[n_param_max-2+jj,ii],5)+' '+ $
                    Strjoin(roundoff(model_err[n_param_max-2+2*jj : n_param_max-2+2*jj+1,ii],5),' ')
        if jj eq 0 then model_str+=' '
      endfor
;     and add to the list of previous components:
      model_str_full= (ii gt 0) ? model_str_full+' | '+model_str : model_str
    endfor  ; for ii=0,max_cmp

    return, model_str_full
END

PRO UPDATE_FILE, $
    stokes_q_arr_orig, stokes_u_arr_orig, filen=filen, new_file=new_file, $
    flag_status_inherited=flag_status_inherited, metrics_arr=metrics
    
;   Write the results from the model fitting to a file.

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common flag_status, flag_status  ; stores error code
 
    common indices, cmp_index

    if keyword_set(flag_status_inherited) then begin
;     Save the flag status of this model, which was inherited from its parent model.
      openw,10,filen,append= ~keyword_set(new_file)  
      tmp_str= Strjoin(Replicate(' 0', 2*n_param_max+2))
      Echo,unit=10,strjoin(roundoff(model_type_arr),' ')+' | '+ $
        fill_str(roundoff(0,2),10)+' '+fill_str(roundoff(0,2),10)+' '+ $  ; loglr, AIC
        fill_str(roundoff(0,2),10)+' '+fill_str(roundoff(99,2),10)+' '+ $  ; BIC, chi2_red
        fill_str(roundoff(100,2),10)+' '+fill_str(roundoff(-1),10)+' '+ $  ; chi2_red_max, det_significance
        fill_str(roundoff(flag_status),10)+ $  ; flag_status     
;       Fill ML estimators for parameters (and their errors) with zeroes:
        Strjoin(Replicate(' |'+tmp_str,cmp_index+1)), /nospawn  
      close,10
      return
    endif

;   Calculate the total number of free parameters in all model components: 
    n_param_tot=0
    for ii=1,8 do begin  ; loop over all possible models
      sel=where(model_type_arr eq ii, n_sel)
      if n_sel gt 0 then begin
        Lookup_model,ii, n_param_used=n_param_used
        n_param_tot += n_sel*n_param_used
      endif
    endfor
    metrics= Calc_metrics(n_param_tot,stokes_q_arr_orig,stokes_u_arr_orig)  
;   returns loglr, AIC, BIC, chi2_red, chi2_red_max, det_significance, flag_status

    model_str_full= Find_model_str()

;   Open a filestream that you will use to write information on the fitted data:
    openw,10,filen,append= ~keyword_set(new_file)  

;   Store information on this fitted model:
    echo,unit=10,strjoin(roundoff(model_type_arr),' ')+' | '+ $
;    print,strjoin(roundoff(model_type_arr),' ')+' | '+ $
      fill_str(roundoff(metrics[0],2),10)+' '+fill_str(roundoff(metrics[1],2),10)+' '+ $  ; loglr, AIC
      fill_str(roundoff(metrics[2],2),10)+' '+fill_str(roundoff(metrics[3],2),10)+' '+ $  ; BIC, chi2_red
      fill_str(roundoff(metrics[4],2),10)+' '+fill_str(roundoff(metrics[5]),10)+' '+ $  ; chi2_red_max, det_significance
      fill_str(roundoff(flag_status),10)+' | '+ $  ; flag_status
;      fill_str('0',10)+' | '+ $
      model_str_full, /nospawn

    close,10
END

PRO INITIALIZE_ARRAYS, $
    stokes_q_arr_orig, stokes_u_arr_orig, filen=filen, parent_index=parent_index

;   Fill the arrays model_param, model_err, stokes_q_mod, and stokes_u_mod.

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

;   Initalize all arrays:
    model_param*=0 & model_err*=0
    stokes_q_mod=stokes_q_arr_orig*0 & stokes_u_mod=stokes_u_arr_orig*0
    if parent_index ne -1 then begin
;     Read in the parameters for the parent model, store these in
;     model_param and model_err, and model the Stokes Q and U
;     frequency spectra based on the parent model:
      filestr_arr=Read_file(filen)
      parent_str=filestr_arr[parent_index]
      tmp=Strsplit(parent_str,'|',/extract,count=n_tmp)
      model_type_arr_parent= fix(Strsplit(tmp[0],' ',/extract))
      tmp=tmp[2:n_tmp-1]  ; remove the model_type_arr and metrics information
      max_cmp_parent=n_tmp-3  ; an index, similar to max_cmp
      p_big_parent=dblarr((max_cmp_parent+1)*n_param_max)  ; p_big for the parent model
      for ii=0,max_cmp_parent do begin
        tmp2=Strsplit(tmp[ii],' ',/extract)
        model_param[0:n_param_max-3,ii]=tmp2[2*indgen(n_param_max-2)]
        model_err[0:n_param_max-3,ii]  =tmp2[2*indgen(n_param_max-2)+1]
;       Treat Q_ref and U_ref separately, because they have both a
;       lower and an upper error (other parameters have only one error):
        model_param[n_param_max-2 : n_param_max-1,ii]= tmp2[2*(n_param_max-2)+[0,3]]
        model_err[n_param_max-2 : n_param_max+1,ii]  = tmp2[2*(n_param_max-2)+[1,2,4,5]]
;
        p_big_parent[ii*n_param_max : (ii+1)*n_param_max-1] = model_param[*,ii]
      endfor  ; for ii=0,max_cmp_parent    
;     Calculate stokes_q_mod and stokes_u_mod:
      fnc_name_arr_parent=strarr(max_cmp_parent+1)
      Fill_fnc_name_arr, model_type=model_type_arr_parent, fnc_name=fnc_name_arr_parent
      stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr_parent, p_big_parent, n_param_max, max_cmp=max_cmp_parent)
      stokes_q_mod=stokes_q_u_mod[0:n_chan-1] & stokes_u_mod=stokes_q_u_mod[n_chan:2*n_chan-1]
    endif  ; if parent_index ne -1 
END

PRO SAVE_RESTORE_DATA_MODELS, save=save, restore=restore, infix=infix, mc_sim_path=mc_sim_path
;   Unfortuntaly variables saved in a .sav file are not automatically
;   restored to a common block. Therefore write this program to save
;   or load variables the way you intend them to be handled.

    common input_model, model_type_arr_input, fnc_name_arr_input, p_big_input1, p_big_input2

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    if n_elements(infix) eq 0 then infix=''
    filen=mc_sim_path+'/fs_mock_data'+infix+'.sav'

    if keyword_set(save) then begin
      stokes_q_arr_cb=stokes_q_arr & stokes_u_arr_cb=stokes_u_arr  ; '_cb' refers to use in a Common Block
      vec_p_freq_weighted_cb=vec_p_freq_weighted
      freq_arr_cb=freq_arr & freq_ref_cb=freq_ref 
      freq_step_cb=freq_step & n_chan_cb=n_chan & lambda2_arr_cb=lambda2_arr
      noise_q_arr_cb=noise_q_arr & noise_u_arr_cb=noise_u_arr
      fwhm_rmsf_cb=fwhm_rmsf
;
      model_type_arr_input_cb=model_type_arr_input 
      fnc_name_arr_input_cb=fnc_name_arr_input
      p_big_input1_cb=p_big_input1
      p_big_input2_cb=p_big_input2

      Save,filen=filen,stokes_q_arr_cb,stokes_u_arr_cb, vec_p_freq_weighted_cb, $
        freq_arr_cb, freq_ref_cb, freq_step_cb, n_chan_cb, lambda2_arr_cb, $
        noise_q_arr_cb, noise_u_arr_cb, fwhm_rmsf_cb, $
        model_type_arr_input_cb, fnc_name_arr_input_cb, $
        p_big_input1_cb, p_big_input2_cb
    endif  ; if keyword_set(save)

    if keyword_set(restore) then begin
      Restore,filen=filen

      stokes_q_arr=stokes_q_arr_cb & stokes_u_arr=stokes_u_arr_cb 
      vec_p_freq_weighted=vec_p_freq_weighted_cb
      freq_arr=freq_arr_cb & freq_ref=freq_ref_cb 
      freq_step=freq_step_cb & n_chan=n_chan_cb & lambda2_arr=lambda2_arr_cb
      noise_q_arr=noise_q_arr_cb & noise_u_arr=noise_u_arr_cb
      fwhm_rmsf=fwhm_rmsf_cb
;
      model_type_arr_input=model_type_arr_input_cb  
      fnc_name_arr_input=fnc_name_arr_input_cb
      p_big_input1=p_big_input1_cb
      p_big_input2=p_big_input2_cb
    endif  ; if keyword_set(restore)
END

PRO GENERATE_MOCK_DATA, seed, infix=infix, mc_sim_path=mc_sim_path, no_noise=no_noise

    common input_model, model_type_arr_input, fnc_name_arr_input, p_big_input1, p_big_input2

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    sel_tmp=where(model_type_arr_input ne 0, n_cmp_inj_max)
    p_big = (n_cmp_inj_max eq 2) ? [p_big_input1,p_big_input2] : p_big_input1
    stokes_q_u_mod=Find_stokes_q_u_mod(fnc_name_arr_input, p_big, n_param_max, max_cmp=n_cmp_inj_max-1)
    stokes_q_arr=stokes_q_u_mod[0 : n_chan-1] & stokes_u_arr=stokes_q_u_mod[n_chan : 2*n_chan-1]

    noise_q_arr = freq_arr*0+1D  ; contains the standard deviation of the noise in Stokes Q
    noise_u_arr = noise_q_arr  ; for the mock observations, assume
                               ; that the noise variances in Q and U
                               ; are the same in each frequency channel.

    if ~keyword_set(no_noise) then begin
;     By default, add noise to the simulated data:
      stokes_q_arr += noise_q_arr*Randomn(seed,n_chan)
      stokes_u_arr += noise_u_arr*Randomn(seed,n_chan)
    endif

    vec_p_freq_weighted=dcomplex(stokes_q_arr/noise_q_arr^2,stokes_u_arr/noise_u_arr^2) 

    Save_restore_data_models, /save, infix=infix, mc_sim_path=mc_sim_path
END

PRO FITIT, $  ;'
    meas_arr=meas_arr, cmp_arr=cmp_arr, n_cmp_max=n_cmp_max, $
    src=src, silent=silent, keep99open=keep99open, $  
;   The following parameters are important if you want to generate mock data:
    test=test, mc_index=mc_index, seed=seed, skip_this_fit=skip_this_fit, $
    l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common indices, cmp_index
;   'cmp_index' is the index of the component that you want to fit

    common param_noise_estimate, delta_logl_onesigma

    common flag_status, flag_status  ; stores error code

    common machine_precision, eps_machine

;   Check that either a data set or a test ID have been supplied:
    if n_elements(meas_arr) eq 0 and n_elements(test) eq 0 then begin
      print,' Please supply a data set or the ID for the mock observation'
      print,''
      print,' For example: fitit, meas_arr=meas_arr, cmp_arr=1, n_cmp_max=5'
      print,''
      print,' meas_arr is an array that measures 7 x n_chan columns x rows'
      print,' (n_chan = the number of observed frequency channels)'
      print,' The seven columns in meas_arr are:'
      print,'  (0) observing frequency (in MHz) '
      print,'  (1,2) measured Stokes I flux density and its 1-sigma uncertainty (both in mJy)'
      print,'  (3,4) measured Stokes Q flux density and its 1-sigma uncertainty (both in mJy)'
      print,'  (5,6) measured Stokes U flux density and its 1-sigma uncertainty (both in mJy)'
      print,''
      print,' cmp_arr is an integer array that lists which model components should be fitted.'
      print,'  numbers range from 1-8, see the explanation at the start of fs.pro.'
      print,''
      print,' n_cmp_max is an integer that indicates the maximum number of components that'
      print,'  should be fitted. Source models with up to and including n_cmp_max components'
      print,'  are fitted to the data.'
      print,''
      stop
    endif

    if n_elements(meas_arr) gt 0 then begin
      Read_measurements, meas_arr, stokes_i_arr=stokes_i_arr, noise_i_arr=noise_i_arr, $
        cmp_arr=cmp_arr,n_cmp_max=n_cmp_max
 
      Lookup_model, n_param_max=n_param_max, param_names_full=param_names_full

      if ~keyword_set(keep99open) then openw,99,'fs.mpfitfun-status.dat'

      if n_elements(src) gt 0 then begin
        filen='fs.'+src+'.dat'  ; write information on the fits to this file
        printf,99,'Source '+src
      endif $
      else filen='fs.dat'
    endif $
    else begin
      if n_elements(test) eq 0 then test=100  ; which test do you want to run by default?
      if n_elements(mc_index) eq 0 then mc_index=0
      if n_elements(seed) eq 0 then seed=0l

      Lookup_test_param, $
        test, test_str=test_str, mc_sim_path=mc_sim_path, cmp_arr=cmp_arr, $
        n_cmp_max=n_cmp_max, n_mc_runs=n_mc_runs, $
;       the following parameters are only required when testing a model with two point sources:  
        l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11

      Lookup_model, n_param_max=n_param_max, param_names_full=param_names_full
 
      infix='.'+test_str+'.index'+roundoff(1+mc_index)

      Generate_mock_data, seed, infix=infix, mc_sim_path=mc_sim_path
;     Note: Call Generate_mock_data after calling Lookup_test_param.
;     Reason for this is that the first function uses Save_restore_data_models, /save
;     to write fwhm_rmsf to file, and this parameter is calculated by Lookup_test_param.

      filen=mc_sim_path+'/fs'+infix+'.dat'

      if keyword_set(skip_this_fit) then goto,skip_to_end_fitit
;     When running batch_fitit with ii_start > 0, first you need to run
;     Generate_mock_data ii_start times to get the correct value for the
;     seed, which is then used by Generate_mock data when you reach
;     ii=ii_start.

      if ~keyword_set(keep99open) then $
        openw,99,mc_sim_path+'/fs.'+test_str+'.mpfitfun-status.dat'
    endelse  ; following if n_elements(meas_arr) ne 0

;   You need the following array to calculate the vector E_j (see equation 11 in DS18):
    vec_p_freq_weighted=dcomplex(stokes_q_arr/noise_q_arr^2,stokes_u_arr/noise_u_arr^2) 
;   stokes_q_arr and stokes_u_arr will be modified by the model
;   fitting. Therefore, make backup copies of the data sets:
    stokes_q_arr_orig=stokes_q_arr & stokes_u_arr_orig=stokes_u_arr

;   List which component types and the maximum number of components you want to fit:
    Fill_model_type_grid, cmp_arr=cmp_arr, n_cmp_max=n_cmp_max, n_models=n_models, $
      model_type_grid=model_type_grid, fnc_name_grid=fnc_name_grid, family_tree_grid=family_tree_grid
;   Store the parameters for each component and the associated errors:
    model_param=dblarr(n_param_max,n_cmp_max)  
;   Store the parameters that describe each model component (one row per component): 
    model_err=dblarr(n_param_max+2,n_cmp_max)
;   The errors in Q_ref and U_ref are described by a lower and an upper number,
;   therefore model_err should have two more columns than model_param.

;   Use the following value for delta log likelihood to calculate
;   the limits of the 68% confidence interval:
    delta_logl_onesigma = 1D /2  
;   originally, you used delta_logl_onesigma = chisqr_cvf(1-0.683D,1)/2D
;   chisqr_cvf(..) = 1, see Mathematica for a distribution with one degree of freedom

;   Use 'machar(/double)' to select when matrix elements are smaller than machine precision:
    eps_machine=Find_eps_machine()
    for ii=0,n_models-1 do begin
;     Select which components make up this model, how many parameters
;     each component has, and the function names (in IDL) of the model
;     components:
      model_type_arr=fix(model_type_grid[*,ii])
;     Find 'cmp_index', which is the index of the model component you want to fit:
      cmp_index=Find_max_cmp()
      fnc_name_arr=reform(fnc_name_grid[*,ii])
      flag_status=family_tree_grid[1,family_tree_grid[0,ii]]  ; flag status of the parent model
      if flag_status ne 0 then begin
        family_tree_grid[1,ii]=flag_status  ; inherit flag status from parent model
        flag_status_inherited=1  
;       This flag indicates that you didn't run Fit_model because of the flag status of the parent model.
        goto,skip_to_update_file
      endif $
      else flag_status_inherited=0
;     If the program arrives at this point it automatically means that flag_status eq 0.    

      Initialize_arrays, $
        stokes_q_arr_orig, stokes_u_arr_orig, filen=filen, parent_index=family_tree_grid[0,ii]

;     Don't refit component 1,2, etc. which were already fitted in the parent model.
;     Instead, read off these parameters from the fit of the parent
;     model, then do a joint fit with the parameters for the new
;     component in the daughter model.
;     Subtract the joint fit (described by parameters in the parent
;     model) from the original measurements:
      stokes_q_arr=stokes_q_arr_orig-stokes_q_mod & stokes_u_arr=stokes_u_arr_orig-stokes_u_mod
;     Fit_model uses these residual spectra to find starting
;     parameters to fit the new model component with MPFITFUN.
      Fit_model, stokes_q_arr_orig, stokes_u_arr_orig

      family_tree_grid[1,ii]=flag_status

      skip_to_update_file:
      Update_file, stokes_q_arr_orig, stokes_u_arr_orig, filen=filen, new_file=(ii eq 0), $
                   flag_status_inherited=flag_status_inherited, metrics_arr=metrics

      if flag_status eq 0 and ~keyword_set(silent) then $
;       Plot the model, using all components:
        Plotit, $
          stokes_q_arr_orig, stokes_u_arr_orig, $
          real_data=n_elements(meas_arr) gt 0, src=src, stokes_i_arr=stokes_i_arr, $
          noise_i_arr=noise_i_arr, chi2_red=metrics[3]
;       I had to move Plotit here, because it depends on the value for
;       the reduced chi squared, which is calculated by Update_file
    endfor  ; for ii=0,n_models-1

    if ~keyword_set(keep99open) then close,99
    skip_to_end_fitit:
END
