!> \file mo_sas_global_variables.f90

!> \brief Global variables used for water quality modelling based on TTDs.

!> \details
!> \note: time variables (e.g, max_age, deni) have a unit of timetep, e.g., 
!>        if the master equation is call at daily, then the unit is day
!>        ..................................monthly (30 days), then the unit is month (x 30 days)

!> \authors Tam Nguyen
!> \date June 2020


module mo_sas_global_variables

  implicit none

  type sas_para
     real*8, dimension(:), allocatable              :: stor_age       ! Residence time distribution
     real*8, dimension(:), allocatable              :: conc_age       ! Concentration distribution
     real*8                                         :: ka             ! k (or a) parameter of the powerlaw (or beta)
     real*8                                         :: b              ! beta parameter (b)
     real*8                                         :: half_life      ! half_life of nitrate [timestep] !change to ln(2)/half_life
  end type sas_para

  type(sas_para), dimension(:), allocatable, public     :: sas_sub

  type sas_output
     integer                              :: median_tt 
     integer                              :: median_rt
     real*8                               :: mean_tt
     real*8                               :: mean_rt
     real*8                               :: denitri_amount
     real*8                               :: subNstore
     real*8, dimension(:), allocatable    :: age_rank_discharge
  end type sas_output
  
  type(sas_output), dimension(:), allocatable, public                  :: sas_out

end module mo_sas_global_variables





