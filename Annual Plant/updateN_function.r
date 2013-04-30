#This is the Watkinson 1980 negative frequency-dependent population growth 
#function incorporating demographic stochasticity used by the script 
#'annualplant_2spp_stoch_par.r.'  

    #2-species annual plant model----------------------------------------------
    updateN = function(r_arg, Nself, N1, a_intra, a1){
      newN = (r_arg * Nself) /(1 + a_intra * Nself + a1 * N1 )
      newN=rpois(1,newN)         #demographic stochasticity
      return(newN)
                  }