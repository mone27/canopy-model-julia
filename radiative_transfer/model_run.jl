using DataFrames
using CSV
using Query
using TimeZones
using ProfileView

include("radiative_transfer_model.jl")
include("parameters_hainich.jl")

function read_fluxnet(path)
    df = CSV.File(path)|> 
      DataFrame |>
      @select(:datetime, :SW_IN_F, :SW_DIF, :SW_OUT, :LW_IN_F, :LW_OUT, :TS_F_MDS_1, :TA_F, :zenith) |>
    # just change the name of the columns and in next step give them a corrent value
      @rename(:SW_IN_F => :sw_sky,
              :SW_DIF => :sw_sky_d,
              :LW_IN_F => :lw_sky,
              :TA_F => :t_leaf,       # Temperature of air aproximation for leaf T for now
              :TS_F_MDS_1 => :t_soil, # Temp soil at 2cm depth
              :SW_OUT => :sw_out,
              :LW_OUT => :lw_out   ) |>
      @mutate(sw_sky_b = _.sw_sky - _.sw_sky_d,  # SW_IN_F is the total radiation direct + diffuse
              t_leaf = 273.15 + _.t_leaf,       # Temperature of air aproximation for leaf T for now
              t_soil = 273.15 + _.t_soil ) |> # Temp soil at 2cm depth            
      DataFrame
      df.datetime = ZonedDateTime.(df.datetime, "yyyy-mm-dd HH:MM:SSzzzz")
      
      return df
  end

  df = read_fluxnet("data/fluxnet_hainich_with_zenith.csv")

  @time out = radiative_transfer_over_input(df, params)