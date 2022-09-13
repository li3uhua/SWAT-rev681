!>       \file mo_sas_master_equation.f90

!>       \brief solve master equation

!>       \details solve master equation to get RTD, TTD, concentration in the outflow

!>       \authors Tam Nguyen

!>       \date June 2020

module mo_sas_master_equation

  use parm,                              only : sas_para
  use mo_sas_utils,                      only : cumsum,        &
                                                cdfbeta,       &
                                                eval_sas

  implicit none

contains

  ! ------------------------------------------------------------------

  !    NAME
  !        master_equation

  !    PURPOSE
  !>       \brief caculates solute export using the StorAge Selection (SAS) function

  !>       \details caculates solute export using the StorAge Selection (SAS) function

  !    INTENT(IN)
  !>       \param[in] "real*8, dimension(:) :: L1_inflow"  inflow to the subsurface [mm]
  !>       \param[in] "real*8, dimension(:) :: L1_outflow" outflow from the subsurface [mm]
  !>       \param[in] "real*8, dimension(:) :: L1_inc"   solute concentration in the inflow [mg/L]
  !>       \param[in] "real*8, dimension(:) :: max_age"   allowable maximum age in storage [timestep]
  !>       \param[in] "real*8, dimension(:) :: max_old_fraction"   allowable maximum old water fraction in storage [-]
  !>       \param[in] "real*8, dimension(:) :: sas_function"   sas function (1 = powerlaw, 2 = beta)

  !    INTENT(INOUT), OPTIONAL
  !>       \param[in] "real*8, dimension(:) :: L1_outc"  solute concentration in the outflow [mg/L]
  !>       \param[in] "real*8, dimension(:) :: sas"  transit time distribution, residence time distribution [-]
  !>       \param[in] "real*8, dimension(:) :: denitri_amount"  amount of nitrate removed by denitrification [kg N /ha]
  !    HISTORY
  !>       \authors Tam Nguyen

  !>       \date June 2020

  subroutine master_equation(sas,                & ! sas_para variables inflow,   
                             inflow,             & ! inflow to the subsurface (mm)       
                             in_conc,            & ! concentration in the inflow (mg/L) 
                             outflow,            & ! outflow from the subsurface (mm)
                             out_conc,           & ! concentration in the out (mg/L)   
                             max_age,            & ! maximum allowable age in storage
                             max_old_fraction,   & ! maximum allowable olest water fraction in storage
                             sas_function,       & ! StorAge selection function (1=powerlaw, 2=beta)
                             median_tt,          & ! median TT
                             median_rt,          & ! median residence time
                             mean_rt,            & ! mean residence time  
                             mean_tt,            & ! median transit time
                             denitri_amount,     & !
                             sub_n_stor,         & ! subsurface N storage (kg/ha)
                             age_rank_discharge)   ! age of water when loading is 50 percent


    implicit none

    type(sas_para), intent(inout)    :: sas
    real*8,         intent(in)       :: inflow
    real*8,         intent(in)       :: in_conc
    real*8,         intent(in)       :: outflow
    real*8,         intent(out)      :: out_conc
    integer,        intent(in)       :: max_age
    real*8,         intent(in)       :: max_old_fraction
    integer,        intent(in)       :: sas_function
    integer,        intent(out)      :: median_tt          ! median transit time
    integer,        intent(out)      :: median_rt          ! median residence time
    real*8,         intent(out)      :: mean_rt            ! median residence time
    real*8,         intent(out)      :: mean_tt            ! median residence time
    real*8,         intent(out)      :: denitri_amount     ! denitrification (kg N ha-1)
    real*8,         intent(out)      :: sub_n_stor         ! subsurface N storage
    real*8, dimension(:), allocatable, intent(out) :: age_rank_discharge

    !local variables
    integer                             :: i, j                                              !counter for loop
    integer                             :: tmax                                              !maximum age in storage or discharge
    real*8, dimension(:), allocatable   :: sas_orig, deriv_sas                               !sas (derivative) function evaluate over the range of normalized age-ranked storage 
    real*8                              :: residual_discharge, discharge_age_i, store_age_i  !subsurface denitrification rate
    real*8                              :: ka, b, temp, eps                                  !parameter of the powerlaw or beta function
    real*8                              :: oldest_water_fraction, age_1                      !volume of water with tmax               
    real*8, dimension(:), allocatable   :: norm_age_rank_stor, stor_age, conc_age

    eps = epsilon(1.0)

    !get paramter of the powerlaw or beta function
    if (sas_function == 1) then
      ka = sas%ka
      b  = 0.0
    else if (sas_function == 2) then
      ka = sas%ka
      b  = sas%b
    else
      stop
    end if

    !check maximum age in storage
    tmax = size(sas%stor_age)
    
    !Denitrification amount (regarless of how much N is exported in discharge)
    denitri_amount = 0.01 *                                                                      &
                     sum((/inflow, sas%stor_age(:) - (/real*8 :: 0.0, sas%stor_age(1:(tmax - 1))/)/)*  &
                         (/in_conc, sas%conc_age/) * (1.0 - exp(-sas%half_life)))

    !print*, tmax, sas%stor_age(tmax)
    !youngest water in storage
    allocate(stor_age(1))
    stor_age(:) = eval_sas((/sas%stor_age(1)/sas%stor_age(tmax)/), sas_function, ka, b)
    age_1 = max(0.0, inflow - outflow * stor_age(1))

    !update normalized age-ranked storage [0,1]
    allocate(norm_age_rank_stor(tmax + 1))
    norm_age_rank_stor(:) = (/0 + age_1, sas%stor_age(:) + age_1/)/ (sas%stor_age(tmax) + age_1)

    !evaluate the sas function over the normalized age-ranked storage
    allocate(sas_orig(tmax + 1))
    sas_orig(:) = eval_sas(norm_age_rank_stor, sas_function, ka, b)

    !*************************************************************Master equation
    !solve the master equation using Mehthod of Lines and Forward Euler
    deallocate(stor_age)
    allocate(stor_age(tmax + 1))

    !stor_age(:) = (/0.0, sas%stor_age(:)/) + inflow - outflow * sas(:)
    !this solution of the master equation does not ensure that stor_age(:)is a monotonically increasing function

    !update age-ranked storage with inflow (the amount of storage with age < Ti)
    stor_age(:) = (/real*8 :: 0.0, sas%stor_age(:)/) + inflow

    !calculate age_rank discharge (the amount of discharge with age < Ti)
    allocate(age_rank_discharge(tmax + 1))
    age_rank_discharge(:) = outflow * sas_orig(:)

    !derivative of stor_age and age_rank_discharge (the volume of storage or discharge with age Ti)
    stor_age(:) = stor_age(:) - (/real*8 :: 0.0, stor_age(1:tmax)/)

    ! if there is no outflow, everthing in outflow is zero
    if (outflow .le. 1.0e-6) then
      age_rank_discharge(:) = 0.0
    
      !update sas
      sas_orig(:) = 0.0
    else 
      age_rank_discharge(:) =  age_rank_discharge(:) - (/real*8 :: 0.0, age_rank_discharge(1:tmax)/)

      !initialize residual discharge
      residual_discharge = 0.0
  
      do i = 1, tmax + 1
  
        !takes water of this age if residual_discharge > 0.0
        age_rank_discharge(i) = age_rank_discharge(i) + residual_discharge
  
        !update residual discharge
        residual_discharge = 0.0
  
        if (age_rank_discharge(i) .le. stor_age(i)) then
  
          !update storage
          stor_age(i) = stor_age(i) - age_rank_discharge(i)
  
        else
  
          !remaining discharge that needs to be taken from older ages
          residual_discharge = age_rank_discharge(i) - stor_age(i)
  
          !update age rank discharge
          age_rank_discharge(i) = stor_age(i)
          !print*,age_rank_discharge(i)
  
          !update stor_age
          stor_age(i) = 0.0
  
        end if
  
      end do
  
      !If with oldest discharge there is no water for outlow, take from all ages (volume weighted)
      age_rank_discharge(:) = age_rank_discharge(:) + residual_discharge * stor_age(:)/sum(stor_age) 
      stor_age(:) =  stor_age(:) - residual_discharge * stor_age(:)/sum(stor_age)

      !update sas
      sas_orig(:) = cumsum(age_rank_discharge)/sum(age_rank_discharge)
    end if

    !Subsurface N storage
    sub_n_stor = 0.01 * sum(exp(-sas%half_life) * (/in_conc, sas%conc_age(:)/) * stor_age(:))

    !convert back to cummulative sum of stor_age and age_rank_discharge
    stor_age(:) = cumsum(stor_age)
    age_rank_discharge(:) =  cumsum(age_rank_discharge)

    !*************************************************************end Master Equation

    !update solute concentration in each parcel
    allocate(conc_age(tmax + 1))
    conc_age(:) = (/in_conc, sas%conc_age(:)/) * exp(-sas%half_life)

    !calcuate pQ * dT
    allocate(deriv_sas(tmax + 1))
    deriv_sas(:) = sas_orig(:) - (/real*8 :: 0.0, sas_orig(1:tmax)/)

    !solute concentration in the outflow
    out_conc = sum(conc_age(:) * deriv_sas(:))
    !print*,"concentration ", conc_age(:)
    !print*,"deriv_sas ",deriv_sas(:)
    print*,"conc ", out_conc, outflow

    !initialized output (median transit time and residence time)
    median_tt = 1
    median_rt = 1
    mean_tt = 0.0

     
    !calcualte median TT, RT50
    do i = 1, size(sas_orig)

      if (outflow .le. 1.0e-6) then
        median_tt = 0.0
      else
        if (sas_orig(i) < 0.5) median_tt = i

        if (i == 1) then
          mean_tt = mean_tt + sas_orig(i)
        else
          mean_tt = mean_tt + (sas_orig(i) - sas_orig(i-1)) * i !?
        end if
      end if


      if (stor_age(i)/stor_age(tmax + 1) < 0.5) median_rt = i
    end do

    ! calculate mean RT
    mean_rt = stor_age(1)

    do i = 2, tmax + 1
      temp = (stor_age(i) - stor_age(i-1))
      if (temp .gt. eps) mean_rt = mean_rt + (stor_age(i) - stor_age(i-1)) * i
    end do

    mean_rt = mean_rt/stor_age(tmax + 1)

    !update TTD and CTD
    deallocate(sas%stor_age)
    deallocate(sas%conc_age)
  
    allocate(sas%stor_age(tmax + 1))
    allocate(sas%conc_age(tmax + 1))   
    sas%stor_age(:) = stor_age(:)
    sas%conc_age(:) = conc_age(:)

  end subroutine master_equation

end module mo_sas_master_equation


























