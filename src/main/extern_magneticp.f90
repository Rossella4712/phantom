module extern_magneticp


implicit none
!default values, can be changed in the input file
 real, public :: dipole_moment = 0.          !6.3*(10**36)     ! magnetic dipole moment of the star [ G*cm^3]
 real, public :: dipole_angle = 0.                    ! angle between stellar magnetic dipole axis and disc axis
 real, public :: cos_dipole_angle = 1.
 real, public :: sin_dipole_angle = 0.
 
 real :: dipole_moment_cgs 

 public :: update_magforce_leapfrog
 public :: get_magnetic_force
 public :: write_options_magforce, read_options_magforce
 private
	
 contains 
	
 ! compute the precession frequency `Omegap'
 subroutine calc_omega_p(xyzi, rhoi, spsoundi, accradius1, Omega_p)  
 use physcon, only: pi, gg
		
 real, intent(in)  :: xyzi(3)         ! Particles position
 real, intent(in)  :: rhoi, spsoundi
 real, intent(in)  :: accradius1
 real, intent(out) :: Omega_p(3)       ! Precession frequency

 !real  :: Omega_k, H_R              
 real :: D 
 real :: xsink,ysink,zsink, modr   
 real :: rsinkgas(3)             ! sink-particle distance 
 real :: mu_vers(3)
 
 xsink = 0.0 
 ysink = 0.0
 zsink = 0.0
 		
 rsinkgas(1) = xyzi(1) - xsink
 rsinkgas(2) = xyzi(2) - ysink
 rsinkgas(3) = xyzi(3) - zsink
 		
 modr = sqrt(dot_product(rsinkgas,rsinkgas))  
 	
 mu_vers(1) = sin_dipole_angle
 mu_vers(2) = 0.
 mu_vers(3) = cos_dipole_angle
 		
 !Omega_k = sqrt(gg*M_star / modr**3)    
 !H_R =  spsoundi / Omega_k
 !D = max(sqrt((r/R_in)**2 - 1) , sqrt(2*H_R/R_in))
 D = modr/accradius1
	
 Omega_p(:) = ( dipole_moment**2 / ((pi**2)*(modr**7)*D*spsoundi*rhoi))*cos_dipole_angle*mu_vers(:)
  


 end subroutine calc_omega_p
 	
 	
 subroutine get_magnetic_force(r, vel, rhoi, spsoundi, accradius1, vcrossomega)
 use vectorutils, only:cross_product3D
 
 real, intent(in)  :: r(3),vel(3)
 real, intent(in)  :: rhoi, spsoundi  
 real, intent(in)  :: accradius1    
 real, intent(out) :: vcrossomega(3)
 real :: Omega_p(3)
 
 call calc_omega_p(r, rhoi, spsoundi, accradius1, Omega_p)
 call cross_product3D(vel,Omega_p,vcrossomega)
 end subroutine get_magnetic_force
 
 
 subroutine update_magforce_leapfrog(vhalfx,vhalfy,vhalfz,fxi,fyi,fzi,rhoi,spsoundi,vcrossomega,dt,xi,yi,zi,accradius1)
 use vectorutils, only : cross_product3D,matrixinvert3D
 use io, only : fatal,warning
		

 real, intent(in)    :: vhalfx,vhalfy,vhalfz
 real, intent(inout) :: fxi,fyi,fzi          !sph force
 real, intent(in)    :: rhoi,spsoundi
 real, intent(out)   :: vcrossomega(3)
 real, intent(in)    :: dt,xi,yi,zi
 real, intent(in)    :: accradius1
		
 integer :: ierr
 real :: dton2
 real :: A(3),v1(3),Omega_p(3) 
 real :: Rmat(3,3),Rinv(3,3)
 real :: xyzi(3)
 
 !print *, xi, yi, zi 
 
 xyzi(1) =  xi
 xyzi(2) =  yi
 xyzi(3) =  zi
 		
 dton2 = 0.5*dt
 
 call calc_omega_p(xyzi, rhoi, spsoundi, accradius1, Omega_p)
 A(1)  = vhalfx + dton2*fxi
 A(2)  = vhalfy + dton2*fyi
 A(3)  = vhalfz + dton2*fzi
		
 Rmat = reshape((/1.,              -dton2*Omega_p(3), dton2*Omega_p(2), &
                  dton2*Omega_p(3), 1.,              -dton2*Omega_p(1), &
                  -dton2*Omega_p(2), dton2*Omega_p(1), 1.             /),(/3,3/))
 
                 
 call matrixinvert3D(Rmat,Rinv,ierr)
 if (ierr /= 0) then
 call fatal('extern_magneticp','Error: determinant = 0 in matrix inversion')
 endif
 			
 v1(:) = matmul(A,Rinv)
 call cross_product3D(v1,Omega_p,vcrossomega)
 
 fxi = fxi + vcrossomega(1)
 fyi = fyi + vcrossomega(2)
 fzi = fzi + vcrossomega(3)
end subroutine update_magforce_leapfrog



!write and read options in input file, runtime

subroutine write_options_magforce(iunit)
 use infile_utils, only:write_inopt
 use physcon, only:pi
 !use units,   only:udist,unit_Bfield
 integer, intent(in) :: iunit
 			
 dipole_angle = dipole_angle*(180.0/pi)
 write(iunit,"(/,a)") '# options relating to magnetic precession'
 call write_inopt(dipole_moment_cgs,'dipole_moment','stellar magnetic dipole moment [G*cm^3] ', iunit) 
 call write_inopt(dipole_angle, 'dipole_angle','angle between stellar magnetic dipole axis and disc axis (0 to 180)', iunit)
 dipole_angle = dipole_angle*(pi/180.0)
end subroutine write_options_magforce


subroutine read_options_magforce(name,valstring,imatch,igotall,ierr)
 use io,      only:fatal
 use physcon, only:pi,au
 use units,   only:udist,unit_Bfield
 character(len=*), intent(in)  :: name,valstring
 logical,          intent(out) :: imatch,igotall
 integer,          intent(out) :: ierr
 integer, save :: ngot = 0
 character(len=30), parameter :: label = 'read_options_magforce'
 	
 imatch  = .true.
 igotall = .false.
 		
 select case(trim(name))
 case('dipole_moment')
    read(valstring,*,iostat=ierr) dipole_moment_cgs  
    if (dipole_moment < 0.) then
     call fatal(label,'invalid dipole moment (should be > 0.)')
    else
     
     ! write dipole_moment in code units, considering udist=1., umass=1., utime=1. 
 
     dipole_moment = dipole_moment_cgs*(unit_Bfield/(sqrt(4.*pi)*8.138E+3))*(udist/au)**3
     print *, dipole_moment
     print *, udist
     print *, unit_Bfield
     
    endif
    ngot = ngot + 1
 case('dipole_angle')
    read(valstring,*,iostat=ierr) dipole_angle
    if (dipole_angle > 180. .or. dipole_angle < 0.) then
     call fatal(label,'invalid dipole angle (should be between 0 and 180 degrees)')
    else
     dipole_angle = dipole_angle*(pi/180.0)
     sin_dipole_angle = sin(dipole_angle)
     cos_dipole_angle = cos(dipole_angle)
    endif
 case default 
   imatch = .false.
 end select
   igotall = (ngot >= 1)
end subroutine read_options_magforce

end module extern_magneticp

