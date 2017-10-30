;   Interpret the information in the files fs*.dat (written by Fitit
;   from fs.pro) in an easy-to-understand format. Print this
;   information to the terminal and plot the data together with the
;   best-fitting model (selected by the Bayesian Information
;   Criterion), unless the keyword 'silent' is used.
;
;   Readit can be applied to data from real observations, or to
;   simulated data from mock observations. These two cases require
;   different sets of input parameters.
;   v For analysing real data, specify at least the array meas_arr.
;   v For analysing Monte Carlo simulations, specify at least values for
;     test, mc_index. If you simulate two source that each emit at a
;     single RM, you should specify also values for l2_test11,
;     chi0_2_test11, and rm2_rayleigh_test11.
;   See also the two examples at the end of this preamble.
;
;   Readit skips parent and daughter models if their flag status is
;   not equal to zero, which indicates that something went wrong in
;   the analysis, if the reduced chi squared of the fit is larger
;   than the value stored in chi2_red_max (one of the columns in
;   metrics_arr), or if the fit converged on a spectral index = -6
;   or +3, which means that the fit converged to the edge of the range
;   of allowed values for \alpha.
;
;   INPUT
;   meas_arr: array with seven columns, and a number of rows equal to
;     the number of frequency channels. The seven columns are:
;     (0) frequency (in MHz) 
;     (1,2) Stokes I flux density and its 1-sigma uncertainty (both in mJy)
;     (3,4) Stokes Q flux density and its 1-sigma uncertainty (both in mJy)
;     (5,6) Stokes U flux density and its 1-sigma uncertainty (both in mJy)
;     No defaults.
;     Specify 'meas_arr' if you want to analyse real data. If this
;     array is not specified, then the program assumes you are
;     analysing a Monte Carlo simulation.
;
;   test: specifies which test from a Monte Carlo simulation to analyse
;         (see the list of tests in fs.pro > Lookup_test_param). 
;         Default is test = 100.
;
;   mc_index: specifies which Monte Carlo run from 'test' to analyse.
;             Note: mc_index starts at zero, default = 0.
;
;   chi2_red_max_override: by specifying a value for 'chi2_red_max_override'
;     fits with values for the reduced chi squared < chi2_red_max_override are
;     included in the analysis. If no value is specified for
;     chi2_red_max_override, then only fits with chi2_red_max are
;     selected by fs.pro > Select_reliable_fits.
;
;   If you want to run Readit, but not Batch_readit, on a double point
;   source model, then use test=11 and specify the following parameters:
;   l2_test11: polarized flux density of the second source component
;
;   chi0_2_test11: intrinsic polarization angle \chi_0 of the second
;                  source component.
;
;   rm2_rayleigh_test11: value for the RM of the second source
;                        component, in units of RM_Rayleigh.
;   No defaults.
;
;   src: name of the source being analysed (not required).
;
;   silent: boolean. If set to 1, no information is written to the
;     terminal and no data are plotted. By default = 0.
;
;   summary: stores a string with information reported by
;     Readit. This way, information can be passed back to the program
;     that calls readit.pro > Readit. For more information, see the text 
;     at the beginning of readit.pro.
;     This variable does not have to be specified.
;
;   only_extract_params: specify '/only_extract_params' if you want to read 
;     the metrics data and model parameters but not write the summary file.
;
;   metrics_arr: stores the metrics that are calculated for the model
;     fits, one row per fit. It contains the following information:
;     (0) model index, (1) log likelihood ratio, (2) AIC, (3) BIC, 
;     (4) chi2_red, (5) chi2_red_max, (6) det_significance, (7) flag_status
;     metrics_arr can be used to pass this information on to programs
;     that call Readit.
;
;
;   OUTPUT,EXAMPLE 1:
;   To report on the tenth realisation simulated for test 101 (see fs.pro > Lookup_test_param), use:
;
;    GDL > readit, test=101, mc_index=9
;
;   Readit then displays which models were fitted, together with the
;   fitted parameters of these models. It also plots the data together
;   with the best-fitting model, as selected by the Bayesian
;   Information Criterion.
;
;   The ID that is listed for a model refers to the row number+1 of that
;   model in the list that is created by fs.pro > Fill_model_type_grid.
;   For example, "Model 1" is the first model in this list etc.
;   Then the program explains which model components were fitted,
;   together with the fitted parameters and their 1-sigma uncertainties.
;
;   For each fit, the program lists a number of metrics:
;
;    chi2_red: the reduced chi squared of the fit
;    chi2_red_max: the largest reduced chi squared that would still be
;                  an acceptable fit (see the penultimate paragraph in
;                  Schnitzeler 2018).
;    AIC: the Akaike Information Criterion
;    BIC: the Bayesian Information Criterion
;    det. significance: the significance of the detection, calculated
;         from the log likelihood ratio, and translated to the equivalent
;         detection significance for a 1D Gaussian distribution to
;         make the number easier to interpret.
;    flag_status: the flag status of the fit, as reported by fs.pro > Fitit.
;                 For more information, see the preamble of fs.pro.
;
;   Next, the program lists model-weighted parameters together with
;   their uncertainties, calculated using equations 18 and 19 in
;   Schnitzeler (2018). The progam calculates these model-weighted
;   parameters first using weights based on the AIC. Then it ranks
;   models based on the value of the AIC, and it shows for each model
;   the AIC weight inside round brackets.
;   Then the program repeats this for the BIC.
;   Last, the program ranks models based on the value for the reduced
;   chi squared. In that case the number inside the brackets is the
;   value for the reduced chi squared.
;
;   In some cases the program will show which models were excluded
;   from further analysis, for example, because there were problems
;   with the fit (fit_status not equal to zero), because the fitting
;   procedure converged on a spectral index of -6 or +3 (the most
;   negative and postive values allowed when fitting the spectral
;   index), or because the reduced chi squared of a fit was larger
;   than the maximum allowed reduced chi squared ('chi2_red_max'). 
;   These fits are not used to rank models or to calculate model weights.
;
;   Readit also plots the data together with the model that was
;   preferred by the Bayesian Information Criterion. What is shown in
;   this plot (or these plots) is explained in fs.pro > Plotit.
;
;
;   EXAMPLE 2:
;   Use the following syntax to apply Readit to real data that have
;   been analysed previously with Fitit:
;
;    GDL > readit, meas_arr=meas_arr, src=src, chi2_red_max_override=199
;
;   The measured frequencies, flux densities, and their uncertainties
;   are stored in the array 'meas_arr'. The layout of this array is
;   explained if you type in 'fitit' after the command prompt.
;   'src' tells the program the name of the source, but you
;     don't need to specify this. 
;   'chi2_red_max_override' tells the program to select fits with a
;     reduced chi squared up to 199, instead of chi2_red_max (the default).
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

@./lib/read_file
@./lib/file_line
@./fs

PRO SORT_MODELS, $
    metrics_arr, column_index, use_chi2_red_metric=use_chi2_red_metric, summary=summary, $
    silent=silent, chi2_red_max_override=chi2_red_max_override, meas_arr=meas_arr, src=src, $
    n_models=n_models, infix=infix, mc_sim_path=mc_sim_path

;   'column_index' specifies which criterion (AIC, BIC, or the reduced
;   chi squared) you use to rank models.
;
;   By specifying a value for 'chi2_red_max_override' all fits with
;   values for the reduced chi squared < chi2_red_max_override are
;   selected, instead of only those fits with chi2_red < chi2_red_max.
;  
;   Note: the first column (index 0) of metrics_arr contains the index
;   of each model, starting at zero.

    common observations, stokes_q_arr, stokes_u_arr, vec_p_freq_weighted, freq_arr, freq_ref, $
                         freq_step, n_chan, lambda2_arr, noise_q_arr, noise_u_arr, fwhm_rmsf

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common model_averaging, param_names_full, n_columns, param_arr, param_err_arr, param_used_arr

    common report, report_i_param_err, report_model_param, report_model_err, report_snr, $
                   report_chi2_red, report_ind

    common machine_precision, eps_machine

    case column_index of
      2: metric_str='AIC'
      3: metric_str='BIC'
      4: metric_str='chi2_red'
    endcase

;   Only consider fits to the data that have a reduced chi squared of
;   at most chi2_red_max or chi2_red_max_override, flag_status = 0,
;   and \alpha not equal to -6 or +3:
    Select_reliable_fits, metrics_arr, param_arr, n_param_max, chi2_red_max_override=chi2_red_max_override, $
                          sel=sel, n_sel=n_sel, compl=sel_compl, ncompl=n_sel_compl  ; function in fs.pro
    if n_sel eq 0 then begin
      suffix=' models were selected by Readit > Sort_models..'
      if n_elements(src) gt 0 then $
        print,' '+src+': no '+suffix $
      else $
        print,' No '+suffix            
      print,' Either the reduced chi squared of each fit was too high, all models were flagged, or \alpha was -6 or +3.'
      summary = '-1(1) '+(n_models gt 1 ? strjoin(replicate('0(0) ',n_models-1)) : '')+' | '+ $
                strjoin(replicate('0(0) ',n_param_max))+' || '+ $
                '-1(1) '+(n_models gt 1 ? strjoin(replicate('0(0) ',n_models-1)) : '')+' | '+$
                strjoin(replicate('0(0) ',n_param_max))+' || '+ $
                '-1(1) '+(n_models gt 1 ? strjoin(replicate('0(99) ',n_models-1)) : '')
;     This string has the same structure as when the model fit had
;     no flags and the reduced chi squared of the fit was close to one.
      report_i_param_err=dblarr(8)
      report_i_param_err[6:7]=-1
      report_model_param=model_param*0
      report_model_err=model_err*0
      report_snr=[-1d,-1d]
      return
    endif

;   Store the properties of the best-fitting model:
    ind_sorted=Sort(metrics_arr[column_index,sel])
    best_model= metrics_arr[0,sel[ind_sorted[0]]]  ; index of the best model
    best_model_cmp= model_type_arr[*,sel[ind_sorted[0]]]  ; which components does the best model contain?
    sel_tmp=where(best_model_cmp ne 0, best_model_n_cmp)  ; how many components does the best model contain?
;   Which parameters (and uncertainties) were fitted in the best model?
    tmp_arr= param_arr[*,sel[ind_sorted[0]]]  
    best_model_param_arr= tmp_arr[0 : best_model_n_cmp*n_param_max-1]
    tmp_arr= param_err_arr[*,sel[ind_sorted[0]]]  
    best_model_param_err_arr= tmp_arr[0 : best_model_n_cmp*n_param_max-1]
;   Find the names of the model components:
    best_model_fnc_name_arr= strarr(best_model_n_cmp)
    Fill_fnc_name_arr, model_type=best_model_cmp, fnc_name=best_model_fnc_name_arr
;   and store the reduced chi squared for the best model fit:
    best_model_chi2_red= metrics_arr[4,sel[ind_sorted[0]]]

    model_str=''
    if ~keyword_set(use_chi2_red_metric) then begin
      if ~keyword_set(silent) then print,' Using the '+metric_str+' weights for model averaging:'

;     Calculate the AIC or BIC weights.
;     Note: the array 'weights_arr' has the same number of rows (=models)
;     as metrics_arr, even though not all models are selected by 'sel'.
;     This makes it easy to look up the weight of each model, without
;     having to take into account if the model was selected based on
;     its flag status.
      weights_arr= dblarr(n_elements(metrics_arr[0,*]))
      weights_arr[sel]= exp(-(metrics_arr[column_index,sel]-metrics_arr[column_index,sel[ind_sorted[0]]])/2)
;     Models that are not selected by 'sel' (meaning that flag_status
;     indicated an error) are automatically assigned a weight of zero.
      weights_arr_normalized=weights_arr/total(weights_arr,/double)

;     Only print weights for models selected by 'sel':
      for ii=0,n_sel-1 do $
        model_str+=' '+roundoff(metrics_arr[0,sel[ind_sorted[ii]]])+ $
          '('+roundoff(weights_arr_normalized[sel[ind_sorted[ii]]],3)+')'

;     Calculate model-weighted values for each parameter.
      param_model_weighted=dblarr(n_columns)  ; stores the model-weighted mean of each parameter
      param_err_model_weighted=param_model_weighted  ; and the associated error
      cmp_index=0l  ; index for looping over each of the source paramaters
      for ii=0l,n_columns-1 do begin
        sel2=where(param_used_arr[ii,*], n_sel2) 
        if n_sel2 gt 0 and total(weights_arr[sel2],/double) gt 0 then begin 
;         Note: the parameter you're averaging can have been
;         kept fixed at zero for a number of models under
;         investigation. In model averaging, only include those 
;         models where this is not the case:
          weights_arr_sel2_normalized = weights_arr[sel2]/total(weights_arr[sel2],/double)
;         You don't have to first select only models where
;         flag_status did not point out an error (i.e., models
;         selected with 'sel'): models with such errors are assigned 
;         a weight of zero in weights_arr, see above.

          param_model_weighted[ii]= total(weights_arr_sel2_normalized*param_arr[ii,sel2],/double)

          param_err_model_weighted[ii]= $
            total(weights_arr_sel2_normalized * $
                  sqrt(param_err_arr[ii,sel2]^2 + (param_arr[ii,sel2]-param_model_weighted[ii])^2), $
                  /double)

;         Report the model-averaged values and their uncertainties:
          index_tmp= ii mod n_param_max  ; which type of source component has been fitted?
          if ~keyword_set(silent) then $
            print,' Component ('+roundoff(fix(ii/n_param_max)+1)+') '+param_names_full[index_tmp]+' : '+ $
              roundoff(param_model_weighted[ii],2)+' +/- '+roundoff(param_err_model_weighted[ii],2) 
        endif  ; if n_sel2 gt 0
;       Note: if n_sel2 eq 0, meaning that this parameter was not
;       fitted in any of the model under consideration, then for this
;       parameter the values of both param_model_weighted and
;       param_err_model_weighted are kept at zero. 
;       The model averaged value of a parameter and its uncertainty
;       should also be kept zero if this parameter only occurred in
;       models that carry a weight of zero.
      endfor  ; for ii=0l,n_columns-1   

      if ~keyword_set(silent) and (ii/n_param_max eq fix(ii/n_param_max)) then print,''
    endif $  ; if ~keyword_set(use_chi2_red_metric)
    else $
      for ii=0,n_sel-1 do $
        model_str+=' '+roundoff(metrics_arr[0,sel[ind_sorted[ii]]])+ $
          '('+roundoff(metrics_arr[column_index,sel[ind_sorted[ii]]],3)+')'

;   To ensure that model_str contains the same number of models,
;   independent of the number of models selected by sel and sel_compl, 
;   add models from sel_compl to model_str:
    suffix= ~keyword_set(use_chi2_red_metric) ? '(0)' : '(99)'
;   Giving the models in sel_compl a weight of zero indicates that
;   they carry a probability of zero to explain the data (Wassermann
;   2000), which agrees with why they were not selected.
    if n_sel_compl gt 0 then $
      for ii=0,n_sel_compl-1 do $
        model_str+=' '+roundoff(metrics_arr[0,sel_compl[ii]])+suffix

    if ~keyword_set(silent) then begin
      print,' Model '+roundoff(best_model)+' has the lowest '+metric_str
      print,'   order of the models:'+model_str
      print,''

      if n_sel_compl gt 0 then begin
        print,' The following models were not selected for model weighting: '+ $
                Strjoin(roundoff(metrics_arr[0,sel_compl]),', ') 
        print,''
      endif 
    endif

    summary+=model_str
    if ~keyword_set(use_chi2_red_metric) then begin
      summary+=' |'
      for ii=0,n_columns-1 do $
        summary+=' '+roundoff(param_model_weighted[ii],5)+ $
          '('+roundoff(param_err_model_weighted[ii],5)+')'
      summary+=' || '
    endif

    if metric_str eq 'BIC' and ~keyword_set(silent) then begin
;     Plot the data and model fit for the model with the smallest BIC

      if n_elements(meas_arr) gt 0 then begin
        dimen=Size(meas_arr,/dim)
        if dimen[0] eq 6 then begin
;         In this case the array 'meas_arr' stores freq - I - err_I - Q - err_Q - U - err_U
          stokes_i_arr=reform(meas_arr[1,*])
          noise_i_arr=reform(meas_arr[2,*])
          stokes_q_arr=reform(meas_arr[3,*])
          noise_q_arr=reform(meas_arr[4,*])
          stokes_u_arr=reform(meas_arr[5,*])
          noise_u_arr=reform(meas_arr[6,*])
        endif $
        else begin
;         'meas_arr' stores freq - I - Q - U - sigma_V, used in older programs.
          stokes_i_arr=reform(meas_arr[1,*])
          stokes_q_arr=reform(meas_arr[2,*])
          stokes_u_arr=reform(meas_arr[3,*])
          noise_q_arr=reform(meas_arr[4,*])
          noise_u_arr=noise_q_arr
        endelse

;       Plotit expects model parameters and their uncertainties to be
;       stored in the arrays model_param and model_err,
;       respectively. Therefore, fill these arrays:
        for ii=0,best_model_n_cmp-1 do begin
          model_param[*,ii]=best_model_param_arr[ii*n_param_max : (ii+1)*n_param_max-1]
          model_err[0:n_param_max-1,ii]=best_model_param_err_arr[ii*n_param_max : (ii+1)*n_param_max-1]
;         model_err lists for Q_ref both -Q_err,+Q_err (same for U_ref). Fix this:
          model_err[n_param_max : n_param_max+1,ii]=model_err[n_param_max-1,ii]
          model_err[n_param_max-1,ii]=model_err[n_param_max-2,ii]
          model_err[n_param_max-2,ii]*=(-1)
          model_err[n_param_max,ii]*=(-1)
        endfor
;       Store model_param and model_err in the following arrays, so
;       that they are available also in Readit:
        report_model_param=model_param
        report_model_err=model_err
      endif $  ; if n_elements(meas_arr) gt 0
      else $
;       Read in the noisy mock data, and the parameters of the injected model:
        Save_restore_data_models, /restore, infix=infix, mc_sim_path=mc_sim_path
;       This reads in the arrays 'stokes_q_arr','stokes_u_arr','noise_q_arr', and 'noise_u_arr'.

      eps_machine=Find_eps_machine()  ; required in Plotit
  
      stokes_q_u_mod_tmp= $
        Find_stokes_q_u_mod(best_model_fnc_name_arr, best_model_param_arr, n_param_max, max_cmp=best_model_n_cmp-1)
      stokes_q_mod=stokes_q_u_mod_tmp[0 : n_chan-1] 
      stokes_u_mod=stokes_q_u_mod_tmp[n_chan : 2*n_chan-1]

      erase
      Plotit, $
        stokes_q_arr, stokes_u_arr, $
        real_data=n_elements(meas_arr) gt 0, src=src, stokes_i_arr=stokes_i_arr, $
        noise_i_arr=noise_i_arr, chi2_red=best_model_chi2_red
      report_chi2_red=best_model_chi2_red
      report_ind=ind_sorted[0]
    endif  ; if metric_str eq 'BIC' and ~keyword_set(silent)
END

PRO REPORT_BEST_MODEL, $
    metrics_arr, summary=summary, silent=silent, chi2_red_max_override=chi2_red_max_override, n_models=n_models, $
    meas_arr=meas_arr, src=src, infix=infix, mc_sim_path=mc_sim_path

;   Sort models based on the value of the AIC, BIC, or chi2_red, and calculate
;   Akaike weights and BIC weights for each of the tested models.
    summary=''  ; will store the information calculated by Sort_models

    Sort_models, metrics_arr, 2, summary=summary, silent=silent, $  ; AIC
      chi2_red_max_override=chi2_red_max_override, meas_arr=meas_arr, $
      src=src, n_models=n_models, infix=infix, mc_sim_path=mc_sim_path

    Sort_models, metrics_arr, 3, summary=summary, silent=silent, $  ; BIC
      chi2_red_max_override=chi2_red_max_override, meas_arr=meas_arr, $
      src=src, n_models=n_models, infix=infix, mc_sim_path=mc_sim_path

    Sort_models, metrics_arr, 4, /use_chi2_red_metric, summary=summary, $  ; chi2_red
      silent=silent, chi2_red_max_override=chi2_red_max_override, $
      meas_arr=meas_arr, src=src, n_models=n_models, infix=infix, $
      mc_sim_path=mc_sim_path
END

PRO READIT, $  ;'
    meas_arr=meas_arr, test=test, mc_index=mc_index, chi2_red_max_override=chi2_red_max_override, $
    l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11, $
    src=src, silent=silent, summary=summary, only_extract_params=only_extract_params, $
    metrics_arr=metrics_arr

    common model, model_type_arr, n_param_max, fnc_name_arr, model_param, model_err, stokes_q_mod, stokes_u_mod

    common model_averaging, param_names_full, n_columns, param_arr, param_err_arr, param_used_arr

    common report, report_i_param_err, report_model_param, report_model_err, report_snr, $
                   report_chi2_red, report_ind

;   Check that either a data set or the ID of the test has been specified:
    if n_elements(meas_arr) eq 0 and n_elements(test) eq 0 then begin
      print,' Readit: please specify the data set or the ID for the mock observation'
      print,''
      print,'  for example: readit, meas_arr=meas_arr (for real observations)'
      print,'  or: readit, test=101, mc_index=0 (for Monte Carlo simulations)'  
      print,''
      print,' Typing "fitit" at the command prompt tells you the layout of meas_arr'
      print,' For an overview of which tests have been hard-coded, see fs.pro > Lookup_test_param'
      print,''
      stop
    endif

    if n_elements(meas_arr) gt 0 then begin
;     Apply Readit to real measurements, not simulations

;     Extract parameters like n_chan from the measurements:
      Read_measurements, meas_arr, stokes_i_arr=stokes_i_arr, noise_i_arr=noise_i_arr

      filen = (n_elements(src) gt 0) ? 'fs.'+src+'.dat' : 'fs.dat'  ; read information on the fits from this file
    endif $
    else begin
;     Use Readit to analyse data from Monte Carlo simulations

      if n_elements(test) eq 0 then test=100  ; which test do you want to analyse?
      if n_elements(mc_index) eq 0 then mc_index=0

      Lookup_test_param, $
        test, test_str=test_str, mc_sim_path=mc_sim_path, $
;       The following parameters are only required if you are testing a
;       model with two point sources:
        l2_test11=l2_test11, chi0_2_test11=chi0_2_test11, rm2_rayleigh_test11=rm2_rayleigh_test11

      infix='.'+test_str+'.index'+roundoff(1+mc_index)
      filen=mc_sim_path+'/fs'+infix+'.dat'
    endelse

    Lookup_model, n_param_max=n_param_max, param_names_full=param_names_full
;   param_names_full lists the names of all parameters.

    if ~keyword_set(silent) then print,'--- List of source components from '+filen+' ---'
;   Find out how many models were fitted ('n_models'), and the largest
;   number of parameters in a single component ('n_param_max'):
    filestr_arr=Read_file(filen, n_lines=n_models)
;   Derive the value of n_cmp_max:
    tmp=strsplit(filestr_arr[n_models-1],'|',/extract)
    n_cmp_max=n_elements(tmp)-2

;   (Re)-initialize:
    model_param=dblarr(n_param_max,n_cmp_max) 
    model_err=dblarr(n_param_max+2,n_cmp_max)

;   Store the values of the fitted parameters of each model (one row
;   per model) and their associated uncertainties:
    n_columns=n_param_max*n_cmp_max  ; enough space to store n_param_max values for each source component
    param_arr=dblarr(n_columns,n_models) 
    param_err_arr=param_arr   
;   Note that Q_ref and U_ref have two errors associated with them (lower and upper values); average these two values.
    model_type_arr=intarr(n_cmp_max,n_models)
    metrics_arr=dblarr(8,n_models)  
;   (0) model index, (1) log likelihood ratio, (2) AIC, (3) BIC, 
;   (4) chi2_red, (5) chi2_red_max, (6) det_significance, (7) flag_status
;
;   Keep track which parameter is used in each model:
    param_used_arr=intarr(n_columns,n_models)  ; value=1: parameter is used (default: 0 -> not used)
;   param_used_arr: do not count a parameter if it was fixed to zero
;   in the fitting (see function Fitit in fs.pro). This is not taken into
;   account properly by default: the normalization of w_i (equations
;   18 and 19 in Schnitzeler 2018) also counts the w_i of models where a
;   parameter is fixed to zero, leading to the wrong result.
;   Using param_used_arr to select the correct models avoids this problem.

    for ii=0l, n_models-1 do begin  ; loop over each model in fs.dat
      if ~keyword_set(silent) then print,'Model '+roundoff(ii+1)

      model_str=Strsplit(filestr_arr[ii],'|',/extract)
;     Extract metrics data:
      metrics=double(Strsplit(model_str[1],' ',/extract))  
;     stores log likelihood ratio - AIC - BIC - chi2_red - chi2_red_max - det_significance - flag status
      metrics_arr[*,ii]=[ii+1,metrics]
;     If this daughter model was skipped because its parent model had
;     a non-zero flag_status, the daughter model will have zeroes for
;     the parameter values (and their errors). You can safely skip
;     reading in these parameters:
      if metrics[6] ne 0 and keyword_set(silent) then goto,skip_to_printing_metrics

      jj=0  ; index for looping through the different components in this model
      ctu=1  ; boolean
      type_arr=fix(Strsplit(model_str[0],' ',/extract,count=n_tmp))
      model_type_arr[0:n_tmp-1,ii]=type_arr  ; stores which components were fitted
      while type_arr[jj] ne 0 and ctu eq 1 do begin
;       Find out which type of model was fitted:
        Lookup_model, type_arr[jj], type_str=type_str, param_used=param_used

        params= double(Strsplit(model_str[2+jj],' ',/extract))
;       Prepare three arrays which you will use later for model-averaging parameters:
        param_arr[jj*n_param_max : (jj+1)*n_param_max-1,ii]= $
          [params[2*indgen(n_param_max-1)], params[2*n_param_max-1]]
        param_err_arr[jj*n_param_max : (jj+1)*n_param_max-1,ii]= $
          [params[2*indgen(n_param_max-2)+1], $
;            Qref and Uref have two errors associated with them: average these numbers:
             (-params[2*(n_param_max-2)+1]+params[2*(n_param_max-2)+2])/2d, $
             (-params[2*n_param_max]+params[2*n_param_max+1])/2d] 
        param_used_arr[jj*n_param_max : (jj+1)*n_param_max-1 ,ii]= param_used 

        if ~keyword_set(silent) then begin
;         Summarize information on the parameters that were fitted.
;         Loop over all elements in this source component, including parameters which
;         are kept zero (sigma_RM_external in a Burn slab, for example):
          param_str='RM= '+roundoff(params[0],2)+'('+roundoff(params[1],2)+'), '  ; Initialize for kk=0
          index=2  ; walk through all elements in "param_str".
          for kk=1, n_param_max-1 do begin  
            if param_used[kk] then begin
              param_str+=param_names_full[kk]+'= '+roundoff(params[index],2)
              if kk le n_param_max-3 then begin
                param_str+='('+roundoff(params[index+1],2)+')'
                index+=2
              endif $
              else begin
;               handle Q_ref and U_ref:
                param_str+='('+Strjoin(roundoff(params[index+1 : index+2],2),',')+')'
                index+=3
              endelse
  
              if kk lt n_param_max-1 then param_str+=', '
            endif $  ; if param_used[kk]
            else index+=2
          endfor  ; for kk=1, n_param_max-1
          print,'('+roundoff(jj+1)+') '+type_str+', '+param_str
        endif  ; if ~keyword_set(silent)

        if jj lt n_cmp_max-1 then jj++ else ctu=0 
      endwhile  ; while type_arr[jj] ne 0 and ctu eq 1

      skip_to_printing_metrics:
      if ~keyword_set(silent) then begin
        print,'    chi2_red: '+roundoff(metrics[3],2)+ $
              '    chi2_red_max: '+roundoff(metrics[4],2)+ $
              '  AIC: '+roundoff(metrics[1],2)+ $
              '  BIC: '+roundoff(metrics[2],2)+ $
              '  det. significance: '+roundoff(metrics[5])+ $
              '  flag status: '+roundoff(metrics[6])
        print,''
      endif
    endfor  ; for ii=0l, n_models-1

    if ~keyword_set(only_extract_params) then $
      Report_best_model, metrics_arr, silent=silent, summary=summary, $
        chi2_red_max_override=chi2_red_max_override, $
        n_models=n_models, meas_arr=meas_arr, src=src, mc_sim_path=mc_sim_path, $
        infix=infix
END
