module montecarlo

 use :: gaussian_core, ONLY : potential, virial, is_zero

 implicit none

 contains

  subroutine put_seed(seed_number)
    integer, intent(in) :: seed_number
    integer, dimension(:), allocatable :: seed 
    integer :: i,n 

    call random_seed(size=n)
    allocate(seed(n))
    do i =1,n 
        seed(i) = seed_number
    end do
    call random_seed(put=seed)
 end subroutine put_seed


 subroutine pbc(ri,box)
    real(kind=8),intent(inout) :: ri(:)   ! position of a single particle 
    real(kind=8),intent(in) :: box(:)       ! simulation bix dimension
    integer :: j 
    do j=1,size(box)
        if (ri(j) > box(j)/2) then 
            ri(j) = ri(j) - box(j)
        else if (ri(j) < - box(j)/2) then 
            ri(j) = ri(j) + box(j)
        end if 
    end do 
 end subroutine pbc

 subroutine interaction(position,box,u,w)
    real(kind=8),intent(in) :: position(:,:)
    real(kind=8),intent(in) :: box(:)
    real(kind=8),intent(out) :: u, w 
    real(kind=8) :: rij(size(position,1)), rij_sq
    real(kind=8) :: uij, wij
    logical :: cutoff
    integer :: i,j 

    u = 0.0
    w = 0.0
    do i=1,size(position,2)
        do j=i+1,size(position,2)
            rij = position(:,i) - position(:,j)
            call pbc(rij,box)
            rij_sq = sum(rij**2)
            call is_zero(rij_sq,cutoff)
            if(cutoff) cycle 
            call potential(rij_sq,uij)
            call virial(rij_sq,wij)
            u = u + uij 
            w = w - wij 
        end do 
    end do 
    w = w/size(box)
 end subroutine interaction 

 subroutine oneparticle_energy(ind,ri,position,box,en)
    integer,intent(in) :: ind
    real(kind=8),intent(in) :: ri(:), box(:), position(:,:)
    real(kind=8),intent(out) :: en 
    real(kind=8) :: rij(size(position,1)), rij_sq, uij 
    logical :: cutoff
    integer :: j
    en = 0.0d0 
    do j=1,size(position,2)
      if(j .NE. ind) then 
          rij = ri - position(:,j)
          call pbc(rij,box)
          rij_sq = sum(rij**2)
          call is_zero(rij_sq,cutoff)
          if(cutoff) cycle
          call potential(rij_sq,uij)
          en = en + uij 
      end if 
    end do 
 end subroutine oneparticle_energy

 subroutine deltaE(ind,trial,position,box,de)
    integer,intent(in) :: ind
    real(kind=8),intent(in) :: trial(:), position(:,:), box(:)
    real(kind=8),intent(out) :: de
    real(kind=8) :: particle_en, trail_en

    call oneparticle_energy(ind,position(:,ind),position,box,particle_en)
    call oneparticle_energy(ind,trial,position,box,trail_en)
    de = trail_en - particle_en    
 end subroutine deltaE

 subroutine move(position,count,dr,box,T)
    integer(kind=8),intent(inout) :: count(1)
    real(kind=8),intent(inout) :: position(:,:)
    real(kind=8),intent(in) :: dr, box(:), T 
    real(kind=8) :: delta(size(position,1)), trial(size(position,1)), deltaU
    real :: rn1, rn2(size(position,1)), rn3 
    integer :: ind, i

    call random_number(rn1)
    ind = 1 + int(rn1*size(position,2))
    call random_number(rn2)
    do i=1,size(position,1)
        delta(i) = dr * (2.0*rn2(i) -1.0)
        trial(i) = position(i,ind) + delta(i)
    end do 
    call pbc(trial,box)
    call deltaE(ind,trial,position,box,deltaU)
    call random_number(rn3)
    if(rn3 <= exp(-deltaU/T)) then
        ! print*,ind
        ! print*,position(:,ind)
        position(:,ind) = trial(:)
        ! print*,position(:,ind)
        count(1) = count(1) + 1  
    end if 
 end subroutine move

 subroutine mc_step(position,count,dr,box,T)
    integer(kind=8),intent(inout) :: count(1)
    real(kind=8),intent(inout) :: position(:,:)
    real(kind=8),intent(in) :: dr, box(:), T 
    integer :: nn 

    do nn=1,size(position,2)
        call move(position,count,dr,box,T)
    end do

 end subroutine mc_step

end module montecarlo
