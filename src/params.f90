!> This module contains core parameters and should not need to be changed.
!! Used by steal and linpop for now.
!! @param i64, i32 64 and 32 bit integers
!! @param dp explicit definition of double precision
!! @param NDIM
!! @param NFDIM 
!! @param MAXIND
!! @param MAXFEIND
module params
    integer, parameter :: i64 = selected_int_kind(15)
    integer, parameter :: i32 = selected_int_kind(6)
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: NDIM = 2560 
    integer, parameter :: MAXIND = 45000 
    integer, parameter :: MAXATOM = 26 
    integer, parameter :: NFDIM   = 2*NDIM + 400 
    integer, parameter :: MAXAUTO = 3200 
    integer, parameter :: MAXINDE = MAXAUTO+MAXIND 
    integer, parameter :: NDDIM   = 89 
    integer, parameter :: NPDIM   = 94 
    integer, parameter :: NFLDIM  = 40 
    integer, parameter :: MAXHIST = 4000 
    integer, parameter :: MAXXDAT = 11 
    integer, parameter :: MAXFEIND = 2500 
end module params