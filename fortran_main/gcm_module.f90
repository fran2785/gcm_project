module gaussian_core
 implicit none

 contains

 subroutine potential(rsq,uij)
    real(kind=8),intent(in)  :: rsq
    real(kind=8),intent(out) :: uij
    real(kind=8),parameter :: epsilon=1.0 , sigma=1.0 

    uij = epsilon*exp(-(rsq/sigma))
 end subroutine potential   

 subroutine virial(rsq,wij)
    real(kind=8),intent(in)  :: rsq
    real(kind=8),intent(out) :: wij
    real(kind=8),parameter :: epsilon=1.0 , sigma=1.0
    real(kind=8) :: r, uij 

    call potential(rsq,uij)
    wij = (-2.0*rsq/sigma)*uij     
 end subroutine virial 

 subroutine is_zero(rsq,cut)
    real(kind=8),intent(in) :: rsq 
    real(kind=8) :: rcut
    logical,intent(out) :: cut 

    rcut = 4.0d0
    cut = rsq > rcut**2      
 end subroutine is_zero

end module gaussian_core
