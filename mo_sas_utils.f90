!>       \file mo_sas_utils.f90

!>       \brief some ultilities functions for solving the storage selection function

!>       \details some ultilities functions for solving the storage selection function

!>       \authors Tam Nguyen

!>       \date June 2020

MODULE mo_sas_utils

  implicit none

  public :: cdfbeta

contains

  function cdfbeta (x, a, b)
  ! ------------------------------------------------------------------

  !    NAME
  !        cdfbeta

  !    PURPOSE
  !>       \brief calculate value of the cummulative beta distribution function at x

  !>       \details calculate value of the cummulative beta distribution function at x
  !>        cdfbeta(x,a,b)   = beta(x,a,b)/beta(a,b) = incomplete beta function / complete beta function

  !    INTENT(IN)
  !>       \param[in] "real*8 :: x" input array
  !>       \param[in] "real*8 :: a" parameter of the beta function
  !>       \param[in] "real*8 :: b" parameter of the beta function

  !    RETURN
  !>       \return real*8 :: value of the beta function at x

  !    HISTORY
  !>       \authors Majumder, K.L., Bhattacharjee, G.P. 

  !    REFERENCE
  !>       Majumder, K.L., Bhattacharjee, G.P. (1973).Algorithm AS 63: The incomplete Beta Integral, 
  !>       Applied Statistics, 22(3), 409-411

  !    MODIFIED
  !>       Tam Nguyen (June, 2020)

  implicit none

  real*8, parameter :: acu = 0.1e-14
  real*8            :: ai
  real*8            :: beta
  real*8            :: cdfbeta
  real*8            :: cx
  logical             :: indx
  integer         :: ns
  real*8            :: a
  real*8            :: pp
  real*8            :: psq
  real*8            :: b
  real*8            :: qq
  real*8            :: rx
  real*8            :: temp
  real*8            :: term
  real*8            :: x
  real*8            :: xx

  beta = log(gamma ( a ) ) + log(gamma ( b ) ) - log(gamma ( a + b ))  
  cdfbeta = x

  !Special cases.
  if ( x == 0.0 .or. x == 1.0 ) then
    return
  end if

  !Change tail if necessary and determine S.
  psq = a + b
  cx = 1.0 - x

  if ( a < psq * x ) then
    xx = cx
    cx = x
    pp = b
    qq = a
    indx = .true.
  else
    xx = x
    pp = a
    qq = b
    indx = .false.
  end if

  term = 1.0
  ai = 1.0
  cdfbeta = 1.0
  ns = int ( qq + cx * psq )

  !Use Soper's reduction formula.
  rx = xx / cx
  temp = qq - ai
  if ( ns == 0 ) then
    rx = xx
  end if

  do

    term = term * temp * rx / ( pp + ai )
    cdfbeta = cdfbeta + term
    temp = abs ( term )

    if ( temp <= acu .and. temp <= acu * cdfbeta ) then

      cdfbeta = cdfbeta * exp ( pp * log ( xx ) &
        + ( qq - 1.0 ) * log ( cx ) - beta ) / pp

      if ( indx ) then
        cdfbeta = 1.0 - cdfbeta
      end if

      exit

    end if

    ai = ai + 1.0
    ns = ns - 1

    if ( 0 <= ns ) then
      temp = qq - ai
      if ( ns == 0 ) then
        rx = xx
      end if
    else
      temp = psq
      psq = psq + 1.0
    end if

  end do

  return

  end function cdfbeta
  ! ------------------------------------------------------------------

  !    NAME
  !        cumsum

  !    PURPOSE
  !>       \brief calculate cummulative sum of an array

  !>       \details calculate cummulative sum of an array
  !>       \cumsum[x]   = /x1, x1 + x2,..., x1 + ... + xn/

  !    INTENT(IN)
  !>       \param[in] "real*8 :: x" input array

  !    RETURN
  !>       \return real*8 :: array of cummulative summation of x

  !    HISTORY
  !>       \authors Tam Nguyen

  !>       \date June 2020

  function cumsum(x)
    implicit none
    
    real*8, dimension(:), intent(in)       :: x
    real*8, dimension(size(x))             :: cumsum
    integer                              :: i

    !initialize result
    cumsum = 0.0

    !first element of the output array
    cumsum(1) = x(1)

    !calculate cumulative summation
    if(size(x) > 1) then
       do i = 2,size(x)
          cumsum(i)  = cumsum(i-1) + x(i)
       end do
    end if
   
  end function cumsum

  ! ------------------------------------------------------------------

  !    NAME
  !        eval_sas

  !    PURPOSE
  !>       \brief caculate value of the sas function (beta or powerlaw)

  !>       \details calculate the powerlaw or beta distribution with given the funciton parameters
  !>       \eval_sas[x,ka]   = powerlaw(x,ka) = x ** ka
  !>       \eval_sas[x,ka,b] = beta(x,ka,b) = see the cdfbeta function

  !    INTENT(IN)
  !>       \param[in] "real*8 :: x" input array

  !    RETURN
  !>       \return real*8 :: values of the powerlaw or beta at x

  !    HISTORY
  !>       \authors Tam Nguyen

  !>       \date June 2020

  function eval_sas(x, sas_function, ka, b)
    implicit none

    real*8,     dimension(:),     intent(in)    :: x   
    integer,                    intent(in)    :: sas_function  
    real*8,                       intent(in)    :: ka      !parameter k (or a) of the powerlaw (or beta) function
    real*8,                       intent(in)    :: b       !parameter b of the beta function

    !local variables
    real*8, dimension(size(x))                  :: eval_sas
    integer                                   :: i, ifault

    !check type of the SAS function
    if (sas_function == 1) then        !powerlaw
       goto 10
    else if (sas_function == 2) then   !beta function
       go to 20
    else                               !unknown funciton
       stop
    end if

    !cummulative powerlaw distribution function
    !if the SAS function is the powerlaw, then ka is the k parameter
10  continue
    eval_sas = x(:) ** ka  !cummulative SAS function
    go to 30

    !cummulative beta distribution function
    !if SAS function is the beta function, ka and b are a, b
20  continue
    do i = 1, size(x)
       eval_sas(i) = cdfbeta(x(i), ka, b)
    end do 

30  continue 

  end function eval_sas

!***************************
  subroutine add_vector(a, b, temp)

    implicit none
    real*8, dimension(:), intent(inout) :: a, b
    real*8, dimension(:), allocatable, intent(out) :: temp

    integer :: n, m

    n = max(size(a), size(b))
    m = min(size(a), size(b))

    allocate(temp(n))

    temp(:) = 0.0
 
    temp(1:m) = a(1:m) + b(1:m)

    if (size(a) < n) temp((m + 1):n) = a(size(a)) + b((m + 1):n)
    if (size(b) < n) temp((m + 1):n) = b(size(b)) + a((m + 1):n)

  end subroutine 


END MODULE mo_sas_utils






