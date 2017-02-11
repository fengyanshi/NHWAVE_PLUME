MODULE PARAM
# if defined (PARALLEL)
     use mpi
# endif
     ! define precision
# if defined (DOUBLE_PRECISION)
     integer, parameter :: SP=8
# if defined (PARALLEL)
     integer, parameter :: MPI_SP=MPI_DOUBLE_PRECISION
# endif
# else
     integer, parameter :: SP=SELECTED_REAL_KIND(6,30)
# if defined (PARALLEL)
     integer, parameter :: MPI_SP=MPI_REAL
# endif
# endif

     ! define parameters
     real(SP), parameter :: pi=3.141592653
     real(SP), parameter :: Small=1.0e-16
     real(SP), parameter :: Large=10000000.0
     real(SP), parameter :: Grav=9.81
     real(SP), parameter :: Zero=0.0
     real(SP), parameter :: One=1.0
     real(SP), parameter :: Rho0=1000.0
     real(SP), parameter :: RhoA=1.20
     real(SP), parameter :: Kappa=0.41

END MODULE PARAM  
