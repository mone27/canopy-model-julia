#' Entry point for the radiative transfer model

include("longwave.jl")
include("shortwave.jl")
include("calc_parameters.jl")

using DataFrames


#' Radiative transfer model step
#'
#' This the core routine of the radiative trasfer model. It calls all the models function
#'
#' @param input A data frame row (or list) containing at least the following elements
#' - datetime
#' - sw_sky_b direct beam shortwave radiation incoming
#' - sw_sky_d diffuse shortwave radiation incoming
#' - lw_sky longwave radiation incoming
#'
#' @param p a list of the model parameters containing at least the following elements
#' max LAI value in the summer
#' min_LAI min value of LAI during winter, it is an aproximation that consider the total Plant Area Index as LAI
#' leaf_out day leaves start in spring
#' leaf_full day leaves reach max LAI
#' leaf_fall day leaves start to fall
#' leaf_fall_complete day all leaves are fallen
#'
#' lat latidude
#' lon longitude
#'
#' rho_leaf Reflencance of leaf
#' tau_leaf trasmissivity of leaf
#' omega_leaf scattering coefficient of leaf
#' clump_OMEGA canopy clumping coefficient
#' alb_soil_b soil albedo direct beam
#' alb_soil_d soil albedo diffuse
#'
#' em_leaf emittivity of leaves
#' em_soil emittivity of soil
#'
#' @return One row data DataFrame with
#' TODO document here

# The Kd in the Two Stream model has a different value
Kd_2stream = get_two_stream_Kd() # This is a costant value that depends only on the leaf angle distribution
function radiative_transfer_model_step(input, p)

  omega_leaf = p[:rho_leaf] + p[:tau_leaf]  
  
  LAI = get_day_LAI(input[:datetime], p[:max_LAI], p[:min_LAI], p[:leaf_out], p[:leaf_full], p[:leaf_fall], p[:leaf_fall_complete])
    zenith = input[:zenith]

    Kb = get_Kb(zenith)
    Kd = get_Kd(LAI)
    beta = get_beta(p[:rho_leaf], p[:tau_leaf])
    beta0 = get_beta0(zenith, Kb, Kd_2stream, omega_leaf)

    shortwave = shortwave_radiation(input[:sw_sky_b], input[:sw_sky_d], LAI, Kb, Kd_2stream, beta, beta0 , omega_leaf,
                                     p[:clump_OMEGA], p[:alb_soil_b], p[:alb_soil_d])
    longwave = longwave_radiation(input[:lw_sky], LAI, input[:t_leaf], input[:t_soil], Kb, Kd, p[:em_leaf], p[:em_soil])

    LAI_sunlit = get_LAI_sunlit(LAI, Kb, p[:clump_OMEGA])

    # values calculated during model run outputed to give more info about the model
    interm_params = DataFrame(:LAI=>LAI, :LAI_sunlit=>LAI_sunlit, :Kb=>Kb, :Kd=>Kd, :beta=>beta, :beta0 => beta0)

    return [shortwave longwave interm_params]

end




#' Runs radiative transfer model over input data
function radiative_transfer_over_input(input, params) 

# maybe not the most efficient stuff, think to preaccolate a matrix
 vcat([radiative_transfer_model_step(row, params) for row in eachrow(input)]...)

end
