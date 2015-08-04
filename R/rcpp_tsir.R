tsir <- function(model_params, model_options, sim_params, summary, init, data, seed){
  .Call( "rcpp_tsir", PACKAGE = "tsir", model_params, model_options,
        sim_params, summary, init, data, seed);
}

