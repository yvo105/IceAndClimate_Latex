#%%
import datetime as dt
import numpy as np
import SEB_functions as SEBf
import matplotlib.pyplot as plt

station = "AWS14"
dtime = "daily"     #'daily' or 'hourly'

# Options Greenland:
# S5 27-08-2003 - 24-08-2023
# S6 29-08-2003 - 21-08-2023
# S9 23-08-2000 - 21-08-2023
# S10 17-08-2010 - 17-04-2026
# Options Antarctica:
# AWS14 21-01-2009 - 01-01-2023
# AWS15 21-09-2009 - 24-06-2014
# AWS17 19-02-2011 - 10-03-2016
# AWS18 25-11-2014 - 02-12-2022
    
if station == "S5" or station == "S6" or station == "S9" or station == "S10":
    SEBdata = SEBf.SEB_data(FileName="Pangaea-Data/GRL_"+ station +"_"+dtime+"_all.csv")
elif station == "AWS14" or station == "AWS15" or station == "AWS17" or station == "AWS18":     
    SEBdata = SEBf.SEB_data(FileName="Pangaea-Data/ANT_"+ station +"_"+dtime+"_all.csv")


#%%
# Example plots of the SEB
# OPTIONS HERE ARE
# Typical yearly cycle -> AvgMonth
# Monthly averages     -> Monthly
# Daily averages       -> Daily

print("some example plots")

PlotType = "Monthly"

if PlotType == "AvgMonth":
    # Example of plotting the typical yearly cycle
    MyFunc = SEBf.get_avg_monthly_value
    Range  = [0, 12.5]
    Label  = "Month"
    Xdata  = np.arange(14)
elif PlotType == "Monthly":
    # Example of montly averages
    MyFunc = SEBf.get_monthly_average
    Range  = [SEBdata.DateTime[0], SEBdata.DateTime[-1]]
    Label  = "Year"
    # do dummy call to function to get the date-times for the x-axis
    Xdata  = MyFunc(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]
elif PlotType == "Daily":
    # Example of montly averages
    MyFunc = SEBf.get_daily_average
    # restrain ourselves to one year, I take here 2016
    Range  = [dt.datetime.fromisoformat("2012-01-01"), dt.datetime.fromisoformat("2013-01-01")]
    Label  = "Date"
    # do dummy call to function to get the date-times for the x-axis
    Xdata  = MyFunc(SHF, SEBdata.DateTime, GiveDatesBack=True)[1]


MonthSWdownS = MyFunc(SWdown, SEBdata.DateTime, PrintInfo=True)
MonthSWupS   = MyFunc(SWup, SEBdata.DateTime)

# however, it is easier to do the call to the SEBf-function in the plot call...

fig, axs = plt.subplots(2, sharex=True)
fig.suptitle(station)

# upper figure are the radiative fluxes
axs[0].plot(Xdata, MonthSWdownS, 'b', linewidth=0.5, label="$SW_{down}$")
axs[0].plot(Xdata, -MonthSWupS, 'b:', linewidth=0.5, label="$SW_{up}$")

axs[0].plot(Xdata, MyFunc(LWdown, SEBdata.DateTime), 'r', linewidth=0.5, label="$LW_{down}$")
axs[0].plot(Xdata, MyFunc(-LWup, SEBdata.DateTime), 'r:', linewidth=0.5, label="$LW_{up}$")
axs[0].plot(Xdata, MyFunc(SWnet+LWnet, SEBdata.DateTime), 'k', linewidth=0.5, label="$R_{net}$")
axs[0].set_ylabel("Energy flux (W/m2)")
axs[0].legend(loc='lower left') # well, no spot is nice
axs[0].grid(True)

axs[1].plot(Xdata, MyFunc(SWnet, SEBdata.DateTime), 'b', linewidth=0.5, label="$SW_{net}$")
axs[1].plot(Xdata, MyFunc(LWnet, SEBdata.DateTime), 'r', linewidth=0.5, label="$LW_{net}$")
axs[1].plot(Xdata, MyFunc(SHF, SEBdata.DateTime), 'seagreen', linewidth=0.5, label="$SHF$")
axs[1].plot(Xdata, MyFunc(LHF, SEBdata.DateTime), 'orange', linewidth=0.5, label="$LHF$")
axs[1].plot(Xdata, MyFunc(Gs, SEBdata.DateTime), 'grey', linewidth=0.5, label="$Gs$")
axs[1].plot(Xdata, MyFunc(MeltS, SEBdata.DateTime), 'purple', linewidth=0.5, label="$M$")
axs[1].plot(Xdata, MyFunc(Residu, SEBdata.DateTime), 'k', linewidth=0.5, label="Residue SEB model")

axs[1].set_ylabel("Energy flux (W/m2)")
axs[1].legend(loc='lower left') # Again, no spot is nice
axs[1].set_xlabel(Label)
axs[1].set_xlim(Range)
axs[1].grid(True)

#%% ALBEDO SENSITIVITY

### FUNCTIONS ###

def SW_net (alpha, SWdown):
    """
    Calculate the net shortwave radiation absorbed by the surface.
    
    Parameters:
    alpha       : surface albedo (unitless)
    SWdown      : Incoming shortwave radiation (W/m^2).
    
    Returns:
    SWnet       : Net absorbed shortwave radiation (W/m^2).
    """
    return (1 - alpha) * SWdown

def LW_up (Tsurf, sigma=5.67e-8):
    """
    Calculate the upward longwave radiation using Stefan-Boltzmann law.
    
    Parameters:
    Tsurf       : Surface temperature (K).
    sigma       : Stefan Boltzmann constant (W/m^2K^4). 
    
    Returns:
    LWup        : Upward longwave radiation (W/m^2).
    """
    return sigma * Tsurf**4

def SHF(Tsurf, T2m, U10m, cs):
    """
    Calculate the sensible heat flux (SHF).
    
    Parameters:
    Tsurf       : Surface temperature (K).
    T2m         : Temperature at 2 meters (K).
    U10m        : Wind speed at 10 meters (m/s).
    cs          : Exchange coefficient for SHF.
    
    Returns:
    SHF         : Sensible heat flux (W/m^2).
    """
    return cs * U10m * (T2m - Tsurf)

def LHF(QSATsurf, QSAT2m, RH2m, U10m, cl):
    """
    Calculate the latent heat flux (LHF).
    
    Parameters:
    QSATsurf    : Saturated specific humidity at the surface (kg/kg).
    QSAT2m      : Saturated specific humidity at 2m above the surface (kg/kg).
    RH2m        : Relative humidity at 2m above the surface.
    U10m        : Wind speed at 10 meters (m/s).
    cl          : Exchange coefficient for LHF.
    
    Returns:
    LHF         : Latent heat flux (W/m^2).
    """
    
    Lhf = cl * U10m * (RH2m * QSATsurf - QSAT2m)
    
    return  Lhf


def energy_balance_linearized(Tsurf_adj, Tsurf_est, SWnet_adj, LWdown_adj, SHFdown_mod, cs, cl, U10m_mod, Tsurf_mod, LHFdown_mod, Gsup_mod, P, RH2m_mod, T2m_mod, sigma=5.67e-8):
    """
    Linearized energy balance equation as a function of the adjusted surface temperature.
    
    Parameters:
    Tsurf_adj               : Adjusted surface temperature (K).
    Tsurf_est               : Initial estimated surface temperature (K).
    SWnet_adj               : Adjusted net shortwave radiation (W/m^2).
    LWdown_adj              : Adjusted incoming longwave radiation (W/m^2).
    SHFdown_mod             : Sensible heat flux from SEB model (W/m^2).
    cs                      : SHF exchange coefficient.
    cl                      : LHF exchange coefficient.
    U10m_mod                : Wind speed at 10 meters (m/s).
    Tsurf_mod               : Surface temperature from SEB (K).
    LHFdown_mod             : Latent heat flux from SEB model (W/m^2).
    Gsup_mod                : Ground heat flux from SEB model (W/m^2).
    P                       : Air pressure (hPa).
    RH2m_mod                : Relative Humidity at 2m (%).
    T2m_mod                 : Air temperature at 2m (K).
    sigma                   : Stefan Boltzmann constant (W/m^2K^4). 
    
    Returns:
    Msurf_adj               : Melt Energy (W/m^2).
    LWup_adj                : Adjustend outward longwave radiation (W/m^2).
    SHFdown_adj             : Adjusted downward sensible heat flux (W/m^2).
    LHFdown_adj             : Adjusted downward latent heat flux (W/m^2).
    """
    
    if Tsurf_adj > 273.16:  
        # cap the surface temperature at melting threshold
        Tsurf_adj = 273.16  
    
    # upward LW_up based on linearised Stefan-Boltzmann law
    LWup_adj = sigma * Tsurf_est**4 + 4 * sigma * Tsurf_est**3 * (Tsurf_adj - Tsurf_est)
    
    # downward SHF based on linearisation
    SHFdown_adj = SHFdown_mod + cs * U10m_mod * (Tsurf_mod - Tsurf_adj)
    
    
    # downward LHF based on linearisation
    Qsat_2m = ClausiusClapeyron(T2m_mod, P)
    Qsat_surf = ClausiusClapeyron(Tsurf_adj, P)
    LHFdown_adj = LHF(Qsat_surf, Qsat_2m, RH2m_mod, U10m_mod, cl*1e3)
    
    
    # Surface energy budget
    Msurf_adj = SWnet_adj + LWdown_adj - LWup_adj + SHFdown_adj + LHFdown_adj + Gsup_mod
    
    return Msurf_adj, LWup_adj, SHFdown_adj, LHFdown_adj


def find_adjusted_surface_temp(Tsurf_est, SWnet_adj, LWdown_adj, SHFdown_mod, cs, cl, U10m_mod, Tsurf_mod, LHFdown_mod, Gsup_mod, P, RH2m_mod, T2m_mod, sigma=5.67e-8, tol=1e-5, max_iter=100):
    """
    Solve for the adjusted surface temperature using the linearized SEB equation.
    
    Parameters:
    Tsurf_est           : Initial estimate of surface temperature (K).
    SWnet_adj           : Adjusted net shortwave radiation (W/m^2).
    LWdown_adj          : Adjusted incoming longwave radiation (W/m^2).
    SHFdown_mod         : Sensible heat flux from SEB (W/m^2).
    cs                  : Constant to adjust SHF with surface temperature change.
    cl                  : Constant to adjust LHF with surface temperature change.
    U10m_mod            : Wind speed at 10 meters (m/s).
    Tsurf_mod           : Surface temperature from SEB (K).
    LHFdown_mod         : Latent heat flux from SEB model (W/m^2).
    Gsup_mod            : Ground heat flux from SEB model (W/m^2).
    P                   : Air pressure (hPa).
    RH2m_mod            : Relative Humidity at 2m (%).
    T2m_mod             : Air temperature at 2m (K).
    sigma               : Stefan Boltzmann constant (W/m^2K^4).
    tol                 : Tolerance for convergence.
    max_iter            : Maximum iterations.
    
    Returns:
    Tsurf_adj_array     : Adjusted surface temperature (K), np.array.
    Msurf_adj_array     : Melt corresponding to adjustment (W/m^2), np.array.
    SHFdown_adj_array   : SHF corresponding to adjustment (W/m^2), np.array.
    LHFdown_adj_array   : LHF corresponding to adjustment (W/m^2), np.array.
    """
    
    # Initialize result arrays
    Tsurf_adj_array = np.ones(len(SWnet_adj))*np.nan
    Msurf_adj_array = np.ones(len(SWnet_adj))*np.nan
    SHFdown_adj_array = np.ones(len(SWnet_adj))*np.nan
    LHFdown_adj_array = np.ones(len(SWnet_adj))*np.nan
    
    # start initial estimate
    Tsurf_adj = Tsurf_est  
    
    for i in range(len(Tsurf_adj_array)):
        
        # check whether nan is present, if yes skip this timestep.
        if (
            np.isnan(SWnet_adj[i]) or
            np.isnan(LWdown_adj[i]) or
            np.isnan(SHFdown_mod[i]) or
            np.isnan(cs[i]) or
            np.isnan(U10m_mod[i]) or
            np.isnan(Tsurf_mod[i]) or
            np.isnan(LHFdown_mod[i]) or
            np.isnan(Gsup_mod[i]) or
            np.isnan(P[i]) or
            np.isnan(RH2m_mod[i])
            ):
            continue  
        
        for _ in range(max_iter):
            
            Msurf_adj, LWup_adj, SHFdown_adj, LHFdown_adj = energy_balance_linearized(
                Tsurf_adj,
                Tsurf_adj,  # Use last adjusted surface temp as estimate for next
                SWnet_adj[i],
                LWdown_adj[i],
                SHFdown_mod[i],
                cs[i],
                cl[i],
                U10m_mod[i],
                Tsurf_mod[i],
                LHFdown_mod[i],
                Gsup_mod[i],
                P[i],
                RH2m_mod[i],
                T2m_mod[i],
                )
            
            
            # if surface temperature exceeds threshold, stop adjusting
            if Tsurf_adj > 273.16:
                Tsurf_adj = 273.16
                break 
            
            if np.abs(Msurf_adj) < tol:
                break  # converged
            
            # Update Tsurf_est for next iteration
            Tsurf_est = Tsurf_adj
            Tsurf_adj = Tsurf_est + Msurf_adj / (4 * sigma * Tsurf_est**3 + (cs[i]) * U10m_mod[i])
    
        # Store the adjusted temperature and melt energy for this timestep
        Tsurf_adj_array[i] = Tsurf_adj
        Msurf_adj_array[i] = Msurf_adj
        SHFdown_adj_array[i] = SHFdown_adj
        LHFdown_adj_array[i] = LHFdown_adj
    
    return Tsurf_adj_array, Msurf_adj_array, SHFdown_adj_array, LHFdown_adj_array

def SHF_CsCalculation (SHF, U10m, T2m, Tsurf):
    """
    Calculates the SHF exchange coefficient.
    
    Parameters:
    SHF     : Observed sensible heat flux (W/m^2).
    U10m    : Wind speed at 10m (m/s).
    T2m     : Air temperature at 2m above ground (K).
    Tsurf   : Surface temperature (K).
    
    Returns:
    cs      : SHF exchange coefficient.
    """
    cs = SHF / (U10m*(T2m - Tsurf))
    return cs

def LHF_ClCalculation (LHF, U10m, Q2m, Qsurf):
    """
    Calculates the LHF exchange coefficient.
    
    Parameters:
    LHF     : Observed latent heat flux (W/m^2).
    U10m    : Wind speed at 10m (m/s).
    Q2m     : Specific humidity at 2m above ground (kg/kg).
    Qsurf   : Specific humidity at surface (kg/kg).
    
    Returns:
    cl      : LHF exchange coefficient.
    """
    cl = LHF / (U10m*(Q2m - Qsurf))
    return cl

def ClausiusClapeyron (T, P, Celsius=False):
    """
    Calculates the saturated specific humidity with the Clausius Clapeyron relation
    
    Parameters:
    T       : Temperature (K).
    P       : Pressure (hPa).
    
    Returns:
    Qsat    : saturated specific humidity (kg/kg)
    """
    m_v = 18.0153  # molecular mass of water vapour (g/mol)
    m_air = 28.9644  # molecular mass of air (g/mol)
    
    if not Celsius:
        T = T - 273.15
    
    c1 = np.where(T <= 0, 22.587, 17.502)
    c2 = np.where(T <= 0, 273.86, 240.97)
    
      
    Qsat = (6.1121 * m_v) / (P * m_air) * np.exp( (c1*T) / (T + c2) )
    
    return Qsat

def get_monthly_mean(Var, DateTime):
    ''' 
    Calculate the monthly mean over all years
    
    Parameters:
    Var         : Variable to be averaged.
    DateTime    : Datetime object corresponding to variable array.
    
    Returns:
    MonthlyMeans: Monthly averages of variable, np.array.
    
    '''
    
    # initiate monthly dictionary
    MonthlyValues = {month: [] for month in range(1, 13)}
    
    # iterate over datetime and retrieve variable values
    for i, date in enumerate(DateTime):
        month = date.month
        MonthlyValues[month].append(Var[i])
    
    # initiate result array
    MonthlyMeans = np.zeros(12)
    
    # calculate monthly average
    for month in range(1, 13):
        MonthlyMeans[month - 1] = np.nanmean(MonthlyValues[month])
    
    return MonthlyMeans

def ClausiusClapeyron_SatVapourPressure (T, Celsius=False):
    """
    Calculates the saturated water vapour pressure with the Clausius Clapeyron relation.
    
    Parameters:
    T       : Temperature (K).
    
    Returns:
    esat    : saturated water vapour pressure (hPa)
    """
      
    if not Celsius:
        T = T - 273.15
    
    c1 = np.where(T <= 0, 22.587, 17.502)
    c2 = np.where(T <= 0, 273.86, 240.97)
    
      
    esat = 6.1121 * np.exp( (c1*T) / (T + c2) )
    
    return esat

def SatMixingRatio_from_SatVapourPressure(P, e_s, Mw=18.0153, Md=28.9644):
    """
    Calculates the saturated mixing ratio from the saturated water vapour pressure. 
    
    Parameters:
    P       : Pressure (hPa).
    q       : Specific humidity (kg/kg).
    Mw      : Molecular mass of water vapour (g/mol).
    Md      : Molecular mass of dry air (g/mol).
    
    Returns:
    w_s     : Saturated mixing ratio.
    """
    epsilon = Mw / Md
    return epsilon * (e_s) / (P - e_s)

def RelativeHumidity_from_SpecificHumidity(T,P,q, Mw=18.0153, Md=28.9644, Celsius=False):
    """
    Calculates relative humidity from specific humidity. 
    
    Parameters:
    T   : Temperature (K).
    P   : Pressure (hPa).
    q   : Specific humidity (kg/kg).
    Mw  : Molecular mass of water vapour (g/mol).
    Md  : Molecular mass of dry air (g/mol).
    
    Returns:
    rh  : Relative humidity.
    """
    if not Celsius:
        T = T - 273.15
    epsilon = Mw / Md    
    w = q / (1-q)
    e_s = ClausiusClapeyron_SatVapourPressure(T)
    w_s = SatMixingRatio_from_SatVapourPressure(P, e_s)
    rh = (w * (epsilon + w_s)) / ((epsilon + w)*w_s)
    
    return rh

#%% #*## SENSITIVITY RUN ###

### SETUP ###

# result dictionary
AlbedoSensitivityDict = {  
    'Albedo': None,
    'Tsurf_adj': None,
    'Msurf_adj': None, 
    'SHFdown_adj': None,
    'LHFdown_adj': None,
}

for key in AlbedoSensitivityDict:
    AlbedoSensitivityDict[key] = {
        'alpha_03': None,
        'alpha_07': None,
        'alpha_085': None
    }
    
AlbedoSensitivityDict["Albedo"] = {
    'alpha_03': 0.3,
    'alpha_07': 0.7,
    'alpha_085': 0.85,
}
Test_ini = 270  # initial guess for surface temperature (K)

## RETRIEVE DATA ##
Gsup_mod = SEBdata.Extract_Variable("GHFup_mod") # SEB model result, calculated using Tsurf_mod
LHFdown_mod = SEBdata.Extract_Variable("LHFdown_mod") # SEB model result, calculated with MO using Tsurf_mod
LWdown_obs = SEBdata.Extract_Variable("LWd") # from observations
Msurf_mod = SEBdata.Extract_Variable("meltE") # SEB model result, calculated using Tsurf_mod
Mcum_mod = SEBdata.Extract_Variable("cumulative_melt") # SEB model result, calculated using Tsurf_mod
P = SEBdata.Extract_Variable("p") # from observations
Q2m_mod = SEBdata.Extract_Variable("q2m") # SEB model result
RH2m_mod = SEBdata.Extract_Variable("rh2m")/100 # SEB model result
SHFdown_mod = SEBdata.Extract_Variable("SHFdown_mod")  # SEB model result, calculated with MO using Tsurf_mod
SHFdown_obs = SEBdata.Extract_Variable("SHFdown")  # from observations
SWdown_obs = SEBdata.Extract_Variable("SWd") # from observations
SWup_obs   = SEBdata.Extract_Variable("SWu") # from observations
T2m_mod    = SEBdata.Extract_Variable("t2m") + 273.15 # modelled air temperature at 2m (similarity theory)
U10m_mod = SEBdata.Extract_Variable('ff10m')  # modelled wind speed at 10m (similarity theory)
Tsurf_mod = SEBdata.Extract_Variable('Ts_mod') + 273.15  # modelled surface temperature in Kelvin.
Tsurf_obs = SEBdata.Extract_Variable('Ts_obs') + 273.15  # modelled surface temperature in Kelvin.


## VARIABLE CALCULATIONS ##
c_s = SHF_CsCalculation(  # calculate SHF exchange coefficient
    SHFdown_mod,
    U10m_mod, 
    T2m_mod,
    Tsurf_mod, 
    )
c_s = np.where(np.isinf(c_s), np.nan, c_s) # replace infinity values with nan

Qsat = ClausiusClapeyron(Tsurf_obs, P)
c_l = LHF_ClCalculation(
    LHFdown_mod,
    U10m_mod,
    Q2m_mod,
    Qsat
)
c_l = np.where(np.isinf(c_l), np.nan, c_l) # replace infinity values with nan

for alphakey in AlbedoSensitivityDict["Albedo"].keys():
    albedo = AlbedoSensitivityDict['Albedo'][alphakey]

    ## APPLY CHANGES IN ALBEDO 
    SWnet_adj = SW_net(albedo, SWdown_obs)

    LWdown_adj = LWdown_obs  # assume that incoming longwave radiation remains 

    # Adjusted surface temperature
    Tsurf_adj, Msurf_adj, SHFdown_adj, LHFdown_adj = find_adjusted_surface_temp(
        Test_ini,
        SWnet_adj, 
        LWdown_adj,
        SHFdown_mod,
        c_s,
        c_l,
        U10m_mod,
        Tsurf_mod,
        LHFdown_mod,
        Gsup_mod,
        P,
        RH2m_mod,
        T2m_mod,
      )
    Msurf_adj = np.where(Msurf_adj > 1e3, np.nan, Msurf_adj)
    AlbedoSensitivityDict["Msurf_adj"][alphakey] = Msurf_adj
    AlbedoSensitivityDict["SHFdown_adj"][alphakey] = SHFdown_adj
    AlbedoSensitivityDict["LHFdown_adj"][alphakey] = LHFdown_adj
    AlbedoSensitivityDict['Tsurf_adj'][alphakey] = Tsurf_adj 
    

_, SEBDateTimeMonths = SEBf.get_monthly_average(AlbedoSensitivityDict['Msurf_adj']['alpha_07'], SEBdata.DateTime, GiveDatesBack=True)


fig, ax = plt.subplots(tight_layout=True)

for key in AlbedoSensitivityDict["Msurf_adj"].keys():
        ax.plot(
        SEBDateTimeMonths,
        SEBf.get_monthly_average(AlbedoSensitivityDict['Msurf_adj'][key], SEBdata.DateTime),
        label=f'Melt Energy with $\\alpha$ = {AlbedoSensitivityDict['Albedo'][key]:.2f}'
        )

ax.plot(
    SEBDateTimeMonths,
    SEBf.get_monthly_average(Msurf_mod, SEBdata.DateTime),
    'k',
    label='Observed Melt Energy'
    )


ax.legend(loc='upper left')
ax.set(
    ylabel='Melt Energy $(W/m^2)$',
    xlabel='Year',
    xlim=(SEBDateTimeMonths[0], SEBDateTimeMonths[-1]),
    title='Melt Energy at AWS 14'

)

#%% LINE PLOT 

Variable = 'Msurf_adj'
VariableName = 'LHF'
VariableUnit = '(W/m^2)'
ObservationName = Msurf_mod


fig, ax = plt.subplots(tight_layout=True)

for key in AlbedoSensitivityDict[Variable].keys():
        ax.plot(
        SEBDateTimeMonths,
        SEBf.get_monthly_average(AlbedoSensitivityDict[Variable][key], SEBdata.DateTime),
        label=f'{VariableName} with $\\alpha$ = {AlbedoSensitivityDict['Albedo'][key]:.2f}'
        )

ax.plot(
    SEBDateTimeMonths,
    SEBf.get_monthly_average(ObservationName, SEBdata.DateTime),
    'k',
    label=f'Observed {VariableName}'
    )


ax.legend(loc='upper left')
ax.set(
    ylabel=f'{VariableName} ${VariableUnit}$',
    xlabel='Year',
    # xlim=(SEBDateTimeMonths[0], SEBDateTimeMonths[-1]),
    ylim=(0,100),
    title=f'{VariableName} at AWS 14'

)


#%% RUNNING MELT SUM 

fig, ax = plt.subplots(tight_layout=True)

for key in AlbedoSensitivityDict[Variable].keys():
        ax.plot(
        SEBdata.DateTime,
        SEBf.get_running_melt_sum(AlbedoSensitivityDict['Msurf_adj'][key], SEBdata.TimeStep,ResetAtNan=False),
        label=f'Melt with $\\alpha$ = {AlbedoSensitivityDict['Albedo'][key]:.2f}'
        )

ax.plot(
    SEBdata.DateTime,
    SEBf.get_running_melt_sum(Msurf_mod, SEBdata.TimeStep,ResetAtNan=False),
    'k',
    label=f'Observed Melt'
    )


ax.legend(loc='upper left')
ax.set(
    ylabel='Accumulated melt (m w.e.)',
    xlabel='Year',
    ylim=(0,60),
)

#%% ALBEDO PLOT

fig, ax = plt.subplots(tight_layout=True)

SWup_monthly, SEBDateTimeMonths = SEBf.get_monthly_average(SWup_obs, SEBdata.DateTime, GiveDatesBack=True)

SWdown_monthly = SEBf.get_monthly_average(SWdown_obs, SEBdata.DateTime)

ax.plot(
    SEBDateTimeMonths,
    SWup_monthly/SWdown_monthly,
    )

ax.set(
    ylabel='Albedo',
    xlabel='Year',
    xlim=(SEBDateTimeMonths[0], SEBDateTimeMonths[-1]),
    ylim=(0,1),

)
ax.grid(alpha=0.7)

#%% LINE PLOT, 1 YEAR ONLY

Variable = 'Msurf_adj'
VariableName = 'LHF'
VariableUnit = '(W/m^2)'
ObservationName = Msurf_mod
ObsMonthAvg, SEBDateTimeMonths = SEBf.get_monthly_average(ObservationName, SEBdata.DateTime, GiveDatesBack=True)
SEBDateTime1Month = [monthindex for monthindex, month in enumerate(SEBDateTimeMonths) if month.year == 2019]
SEBDateTimeMonths = SEBDateTimeMonths[SEBDateTime1Month]


fig, ax = plt.subplots(tight_layout=True)

for key in AlbedoSensitivityDict[Variable].keys():
        VarMonthAvg = SEBf.get_monthly_average(AlbedoSensitivityDict[Variable][key], SEBdata.DateTime)
        ax.plot(
        SEBDateTimeMonths,
        VarMonthAvg[SEBDateTime1Month],
        label=f'{VariableName} with $\\alpha$ = {AlbedoSensitivityDict['Albedo'][key]:.2f}'
        )

ax.plot(
    SEBDateTimeMonths,
    ObsMonthAvg[SEBDateTime1Month],
    'k',
    label=f'Observed {VariableName}'
    )


ax.legend(loc='upper left')
ax.set(
    ylabel=f'{VariableName} ${VariableUnit}$',
    xlabel='Year',
    ylim=(0,50),
    title=f'{VariableName} at AWS 14'

)


#%% #*### BAR PLOT 

## PLOT SETUP
Var = 'Msurf_adj'
VarOrig = Msurf_mod
VarName = 'LHF'
VarUnit = 'W/m^2'
AlbedoName = 'alpha_085'
AlbedoValue = 0.85
fig, ax = plt.subplots(tight_layout=True)
months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Calculate Monthly Means
VarAdj_MonthlyMean = get_monthly_mean(AlbedoSensitivityDict[Var][AlbedoName], SEBDateTimeMonths)
VarOrig_MonthlyMean = get_monthly_mean(VarOrig, SEBDateTimeMonths)


# plot model under new surface temperature
ax.bar(months, VarAdj_MonthlyMean,
       capsize=5,
       alpha=0.7,
       label=f'{VarName} with $\\alpha$ = {AlbedoValue}',
       color='tab:blue',
       )

# plot original model data
ax.bar(months, VarOrig_MonthlyMean,
       capsize=5,
       alpha=0.5,
       color='k',
       label=f'Original {VarName}')

# PLOT SETTINGS

ax.legend(loc='upper left')
ax.set(
    ylabel=f'{VarName} $({VarUnit})$',
    ylim=(0, 50), 
)
ax.grid(alpha=0.7)




