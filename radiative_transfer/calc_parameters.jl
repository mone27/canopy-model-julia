import Dates

#' Simple model to get current LAI
#'
#' Use a simple model to get the LAI of the different day of the year.
#' It has 4 phases:
#'  - Winter: from leaf_fall_complete to leaf_out -> LAI = 0
#'  - Spring: from leaf_out to leaf_full -> linear growth from 0 to max_LAI
#'  - Summer: from leaf_full to leaf_fall -> LAI = max_LAI
#'  - Fall: from leaf_fall to leaf_fall_complete -> linear decrease from max_LAI to 0
#' The 4 paramenters (leaf_out...) are the day of the year
#'
#' @param time a datetime object
#' @param max_LAI max LAI value in the summer
#' @param min_LAI min value of LAI during winter, it is an aproximation that consider the total Plant Area Index as LAI
#' @param leaf_out day leaves start in spring
#' @param leaf_full day leaves reach max LAI
#' @param leaf_fall day leaves start to fall
#' @param leaf_fall_complete day all leaves are fallen
#' 
#' @return LAI Leaf Area Index value for the day of the year
function get_day_LAI(datetime, max_LAI, min_LAI, leaf_out, leaf_full, leaf_fall, leaf_fall_complete)

    yday = Dates.dayofyear(datetime)
    if (yday < leaf_out)  # before leaves are out LAI is min
        return(min_LAI)
    end
    if (yday >= leaf_out & yday < leaf_full ) 
        ndays =  leaf_full - leaf_out # n days of the transition
        return((max_LAI - min_LAI) * (yday - leaf_out) / ndays + min_LAI)
    end
    if (yday >= leaf_full & yday < leaf_fall ) 
        return(max_LAI)
    end
    if (yday >= leaf_fall & yday < leaf_fall_complete ) 
        ndays =  leaf_fall_complete - leaf_fall # n days of the transition
        return((max_LAI - min_LAI) * (leaf_fall_complete - yday) / ndays + min_LAI)
    end
    if (yday >= leaf_fall_complete)
        return(min_LAI)
    end
end
    
    
    
#'  All the following function assumes a SPHERICAL leaves distribution
#' Chapter 2.2


#' Direct beam extiction coefficient

#' @param zenith in degrees
#' @return Kb
get_Kb =function(zenith, max_Kb=20)
    # Eq. 14.29
    Kb = 0.5/cos(deg2rad(zenith)) # extinction coefficient
    Kb = min(Kb, max_Kb) # Prevent the Kb to become too large at low sun angles.
    # The default value of 20 is from the Bonan matlab code script sp_14_03 line 150
    return(Kb)
end

#' Diffuse beam extiction coefficient
#' @param LAI
#' @return Kd
function get_Kd(LAI)
    G_z =  0.5

    # Eq. 14.33
    td =  0
    for z in 0:pi / 18:pi / 4 # make 9 steps from 0 till π/2
        td = td + exp( - G_z / cos(z) * LAI)*sin(z)*cos(z)*(pi / 18)
    end

    # Eq 14.34
    Kd = -log(2 * td)/LAI
    return(Kd)
end


function get_two_stream_Kd()
    # Eq. 14.31
    ross = 0.01 # should be zero but if is zero it mess up the computations.
    # See Bonan matlab code script sp_14_03 line 130
    phi_1 = 0.5 - 0.633 * ross - 0.333 * (ross)^2
    phi_2 = 0.877 * (1 - 2 * phi_1 )
    # Eq 14.80
    Kd =  1 / (( 1 - phi_1/phi_2 * log((phi_1+phi_2)/phi_1) ) / phi_2)
    return(Kd)
end


#' Fraction of diffuse light scattered backward
#' @param rho_leaf
#' @param tau_leaf
#' 
#' @return beta
function get_beta(rho_leaf, tau_leaf) 
    # Derived from equations 14.81 following the book approximation for sperical distribution
    beta = ( 0.625 * rho_leaf +  0.375 * tau_leaf ) / (rho_leaf + tau_leaf)
    return(beta)

end


#' Fraction of direct light scattered backward
#' @param zenith in degrees
#' @param Kb
#' @param Kd
#' @param omega_leaf
#'
#' @return beta0
function get_beta0(zenith, Kb, Kd, omega_leaf)

    # Eq. 14.31
    ross = 0
    phi_1 = 0.5 - 0.633 * ross - 0.333 * (ross)^2
    phi_2 = 0.877 * (1 - 2 * phi_1 )

    G_mu = 0.5 #mu is cos(Z) but G(Z) for sperical leaves distribution is costant
    mu = cos(deg2rad(zenith))

    # Equation 14.84

    #defining commonly used terms
    mphi_1 = mu * phi_1
    mphi_2 = mu * phi_2

    a_s = (
    (omega_leaf / 2) * (G_mu) / (G_mu + mphi_2) *
    (1 - (mphi_1/(G_mu + mphi_2) * log((G_mu + mphi_1 + mphi_2) / mphi_1)))
    )

    beta_0 =  (((Kb + Kd) / Kb) * a_s ) / omega_leaf
    return(beta_0)
end

function get_LAI_sunlit(LAI, Kb, clump_OMEGA)
    # Eq.14.18 integrated in the same way of Eq. 14.12 (also line in Bonan Matlab code line script sp_14_03 line 167)
    LAI_sunlit =  (1 - exp(- clump_OMEGA * Kb * LAI) )/ Kb
    return(LAI_sunlit)
end
    