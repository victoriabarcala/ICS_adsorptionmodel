# this function assumes that the amount of ICS in the column only modifies K linearly

P_column_k <- function (t, state, pars)   # state is a LONG vector
{
  with (as.list(pars),{ 
    
    # unpack state variables
    Paq       <- state[ (0*N+1) : (1*N) ]    # first N elements: Paq
    Pads_slow <- state[ (1*N+1) : (2*N) ]    # next N elements:  Pads_slow
    
    # === transport rates ===
    # note: zero gradient by default at lower boundaries
    
    # dissolved phosphate  
    
    tran.Paq <-  tran.1D(C = Paq, C.up = Pinflow,     # upper boundary: flux 
                         dx = Grid, VF = porosity,    # grid and volume fraction
                         D = Dispers, v = v_adv[t])   # mixing 
    
    
    
    # === reaction rates ===
    # Kinetic adsorption of phosphate
    
    rate_ads <- alpha * (x*(1 - f) * k_ads * Paq - Pads_slow)
    
    # ====mass balances====
    
    dPaq.dt       <- if(v_adv[t] > 0)((tran.Paq$dC - rate_ads*rho/porosity)/(1 +k_ads*f*x*rho/porosity)) else ((-rate_ads*rho/porosity)/(1 +k_ads*f*x*rho/porosity))
    
    
    dPads_slow.dt <-  rate_ads
    
    # concentration of adsorbed P in equilibrium with solution
    Pads_fast     <-   k_ads *x* f * Paq #[mg/g-solid]
    
    
    # depth-integrated state variables
    TotalPads_slow  <- sum (Pads_slow * Grid$dx*rho/porosity)    # [mg-P*cm/ L water]
    TotalPads_fast  <- sum (Pads_fast * Grid$dx *rho/porosity)   # [mg-P*cm/ L water]
    
    #calculations for mass balances
    Paq_outflux <-   Paq[N] * v_adv[t]        # [mg-P *cm/L water/h]
    Paq_influx  <-   Pinflow * v_adv[t]       # [mg-P *cm/L water/h]
    Pretention  <-   Paq_influx - Paq_outflux # [mg-P *cm/L water/h]
    
    
    return(list(c(dPaq.dt, dPads_slow.dt),   Paq_out = Paq[N],          # the time-derivatives, as a long vector
                
                TotalPads_slow = TotalPads_slow,
                TotalPads_fast = TotalPads_fast,
                
               
                # for creating budgets
                Paq_outflux = Paq_outflux,      # [mg-P *cm/L water/h]
                Paq_influx  = Paq_influx ,      # [mg-P *cm/L water/h]
                Pretention  = Pretention,       # [mg-P *cm/L water/h]
                Padsfast=Pads_fast[N], Padsslow=Pads_slow[N]))   
  })
}