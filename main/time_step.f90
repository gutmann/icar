!> ----------------------------------------------------------------------------
!!  Main time stepping module. 
!!  Calculates a stable time step (dt) and loops over physics calls
!!  Also updates boundaries every time step. 
!!
!!  @author
!!  Ethan Gutmann (gutmann@ucar.edu)
!!
!! ----------------------------------------------------------------------------
module time_step
    use data_structures     ! *_type  types
    use microphysics,               only : mp
    use convection,                 only : convect
    use land_surface,               only : lsm
    use wind,                       only : balance_uvw
    use advection,                  only : advect
    use output,                     only : write_domain
    use planetary_boundary_layer,   only : pbl
    use radiation,                  only : rad
    use boundary_conditions,        only : update_pressure
    use io_routines,                only : io_write
    implicit none
    private
    public :: step

    real, dimension(:,:), allocatable :: lastw, currw, uw, vw !> temporaries used to compute diagnostic w_real field
contains
    
    !>------------------------------------------------------------
    !!  Update the edges of curdata by adding dXdt
    !!
    !!  Apply dXdt to the boundaries of a given data array. 
    !!  For these variables, dXdt is a small array just storing the boundary increments
    !!
    !! @param curdata   data array to be updated
    !! @param dXdt      Array containing increments to be applied to the boundaries. 
    !!                  (nz x max(nx,ny) x 4)
    !!
    !!------------------------------------------------------------
    subroutine boundary_update(curdata,dXdt)
        implicit none
        real,dimension(:,:,:), intent(inout) :: curdata
        real,dimension(:,:,:), intent(in) :: dXdt
        integer::nx,nz,ny,i
        
        nx=size(curdata,1)
        nz=size(curdata,2)
        ny=size(curdata,3)

        do i=1,nz
            curdata(1,i,2:ny-1)  = curdata(1,i,2:ny-1)  + dXdt(i,2:ny-1,1)
            curdata(nx,i,2:ny-1) = curdata(nx,i,2:ny-1) + dXdt(i,2:ny-1,2)
            curdata(:,i,1)  = curdata(:,i,1)  + dXdt(i,1:nx,3)
            curdata(:,i,ny) = curdata(:,i,ny) + dXdt(i,1:nx,4)
        enddo
        ! correct possible rounding errors, primarily an issue of clouds...
        where(curdata(1,:,:)<0) curdata(1,:,:)=0
        where(curdata(nx,:,:)<0) curdata(nx,:,:)=0
        where(curdata(:,:,1)<0) curdata(:,:,1)=0
        where(curdata(:,:,ny)<0) curdata(:,:,ny)=0
        
    end subroutine boundary_update
    
    !>------------------------------------------------------------
    !!  Updated all fields in domain using the respective dXdt variables
    !!
    !!  
    !!
    !! @param domain    full domain data structure to be updated
    !! @param bc        boundary conditions (containing dXdt increments)
    !! @param options   options structure specifies which variables need to be updated
    !!
    !!------------------------------------------------------------
    subroutine forcing_update(domain,bc,options)
        implicit none
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        integer::j,ny
        ny=size(domain%p,3)
        !$omp parallel firstprivate(ny) &
        !$omp private(j) &
        !$omp shared(bc,domain)
        !$omp do schedule(static)
        do j=1,ny
            domain%u(:,:,j)=domain%u(:,:,j)+bc%du_dt(:,:,j)
            domain%v(:,:,j)=domain%v(:,:,j)+bc%dv_dt(:,:,j)
            if (.not.options%ideal) then
                domain%p(:,:,j)=domain%p(:,:,j)+bc%dp_dt(:,:,j)
                domain%pii(:,:,j)=(domain%p(:,:,j)/100000.0)**(Rd/cp)
                domain%rho(:,:,j)=domain%p(:,:,j)/(Rd*domain%th(:,:,j)*domain%pii(:,:,j)) ! kg/m^3
            endif
            
            ! these only get updated if we are using the fluxes derived from the forcing model
            if (options%physics%landsurface==kLSM_BASIC) then
                domain%sensible_heat(:,j)  = domain%sensible_heat(:,j)+ bc%dsh_dt(:,j)
                domain%latent_heat(:,j)    = domain%latent_heat(:,j)  + bc%dlh_dt(:,j)
            endif
            
            ! these only get updated if we are using the fluxes (and PBL height) derived from the forcing model
            if (options%physics%boundarylayer==kPBL_BASIC) then
                domain%pbl_height(:,j) = domain%pbl_height(:,j)+ bc%dpblh_dt(:,j)
            endif
            
            ! these only get updated if we are using the fluxes derived from the forcing model
            if (options%physics%radiation==kRA_BASIC) then
                domain%swdown(:,j)  = domain%swdown(:,j)  + bc%dsw_dt(:,j)
                domain%lwdown(:,j)  = domain%lwdown(:,j)  + bc%dlw_dt(:,j)
            endif
            domain%sst(:,j)  = domain%sst(:,j)  + bc%dsst_dt(:,j)
        enddo
        !$omp end do
        !$omp end parallel
        ! v has one more y element than others
        ny=ny+1
        domain%v(:,:,ny)=domain%v(:,:,ny)+bc%dv_dt(:,:,ny)
        ! dXdt for qv,qc,th are only applied to the boundarys
        if (.not.options%ideal) then
            call boundary_update(domain%th,bc%dth_dt)
            call boundary_update(domain%qv,bc%dqv_dt)
            call boundary_update(domain%cloud,bc%dqc_dt)
        endif
        
        ! because density changes with each time step, u/v/w have to be rebalanced as well. 
        ! could avoid this by assuming density doesn't change... but would need to keep an "old" density around
        ! also, convection can modify u and v so it needs to rebalanced
        call balance_uvw(domain,options)
        
        write(*,*) "before diagnostic update in forcing_update"
        write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
        write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
        write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
        write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
        write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
        write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
        write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
        write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
        write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
        write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
        write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
        write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
        write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
        write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
        write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
        write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
        write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
        write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
        write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
        write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
        write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
        write(*,*) "call diagnostic_update in forcing update"
        !call io_write("bc%dqv_dt_bdu.nc","data",bc%dqv_dt)
        call diagnostic_update(domain,options)
        write(*,*) "after diagnostic update in forcing update"
        !call io_write("bc%dqv_dt_adu.nc","data",bc%dqv_dt)
        write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
        write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
        write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
        write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
        write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
        write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
        write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
        write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
        write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
        write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
        write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
        write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
        write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
        write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
        write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
        write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
        write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
        write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
        write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
        write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
        write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
        write(*,*) "end forcing update"
    end subroutine forcing_update
    
    !>------------------------------------------------------------
    !! Update model diagnostic fields
    !!
    !! Calculates most model diagnostic fields such as Psfc, 10m height winds and ustar
    !!
    !! @param domain    Model domain data structure to be updated
    !! @param options   Model options (not used at present)
    !!
    !!------------------------------------------------------------
    subroutine diagnostic_update(domain,options)
        implicit none
        type(domain_type), intent(inout) :: domain
        type(options_type), intent(in) :: options
        
        real :: testlog1, testlog2
        !real, dimension(:,:,:), allocatable :: testlog1, testlog2
        integer :: nx,ny,nz, y, z
        integer :: i1,i2
        
        nx=size(domain%p,1)
        nz=size(domain%p,2)
        ny=size(domain%p,3)
        
        ! update p_inter, psfc, ptop, Um, Vm, mut
        domain%Um = 0.5*(domain%u(1:nx-1,:,:)+domain%u(2:nx,:,:))
        domain%Vm = 0.5*(domain%v(:,:,1:ny-1)+domain%v(:,:,2:ny))
        domain%t  = domain%th * domain%pii
        
        domain%p_inter=domain%p
        call update_pressure(domain%p_inter, domain%z, domain%z_inter, domain%t)
        domain%psfc = domain%p_inter(:,1,:)
        ! technically this isn't correct, we should be using update_pressure or similar to solve this
        domain%ptop = 2*domain%p(:,nz,:) - domain%p(:,nz-1,:)
        
        ! dry mass in the gridcell is equivalent to the difference in pressure from top to bottom
        domain%mut(:,1:nz-1,:) = domain%p_inter(:,1:nz-1,:) - domain%p_inter(:,2:nz,:)
        domain%mut(:,nz,:) = domain%p_inter(:,nz,:) - domain%ptop
        
        if (.not.allocated(lastw)) then
            allocate(lastw(nx-2,ny-2))
            allocate(currw(nx-2,ny-2))
            allocate(uw(nx-1,ny-2))
            allocate(vw(nx-2,ny-1))
            !allocate(testlog1(3000,nx-2,ny-1))
            !allocate(testlog2(3000,nx-2,ny-1))
        endif
        
        ! ----- start surface layer calculations usually done by surface layer scheme ----- !
        write(*,*) "start surface layer calculations"
        if (options%physics%boundarylayer==kPBL_SIMPLE) then
            write(*,*) "calculate surface layer based on log wind profile"
            ! temporary constant
            ! use log-law of the wall to convert from first model level to surface
            currw = karman / log((domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) &
                    / domain%znt(2:nx-1,2:ny-1))
            ! use log-law of the wall to convert from surface to 10m height
            lastw = log(10.0 / domain%znt(2:nx-1,2:ny-1)) / karman
            domain%ustar(2:nx-1,2:ny-1) = domain%Um(2:nx-1,1,2:ny-1) * currw
            domain%u10(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) * lastw
            domain%ustar(2:nx-1,2:ny-1) = domain%Vm(2:nx-1,1,2:ny-1) * currw
            domain%v10(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) * lastw
        
            ! now calculate master ustar based on U and V combined in quadrature
            domain%ustar(2:nx-1,2:ny-1) = sqrt(domain%Um(2:nx-1,1,2:ny-1)**2 & 
                                        + domain%Vm(2:nx-1,1,2:ny-1)**2) * currw
            ! counter is just a variable helping me to detect how much rounds
            ! this subroutine went through
            write(*,*) "Counter: ", counter
            counter = counter + 1
            write(*,*) "Counter: ", counter
        elseif (options%physics%boundarylayer==kPBL_YSU) then
            ! start surface layer calculations introduced by Patrik Bohlinger
            write(*,*) "calculate surface layer based on monin-obukhov similarity theory"

            ! ----- start temporary solution ----- !
            ! use log-law of the wall to convert from first model level to surface
            currw = karman / log((domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) &
                   / domain%znt(2:nx-1,2:ny-1))
            ! use log-law of the wall to convert from surface to 10m height
            lastw = log(10.0 / domain%znt(2:nx-1,2:ny-1)) / karman
            ! calculate ustar = horizontal wind speed scale
            if (counter==1) then
                domain%ustar_new(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1)
            endif
            ! preventing ustar from being smaller than 0.1 as it could be under
            ! very stable conditions, Jiminez et al. 2012
            where(domain%ustar_new(2:nx-1,2:ny-1) < 0.1)
                domain%ustar_new(2:nx-1,2:ny-1) = 0.1
            endwhere
            !domain%ustar(2:nx-1,2:ny-1) = domain%Um(2:nx-1,1,2:ny-1) * currw
            domain%u10(2:nx-1,2:ny-1) = domain%ustar_new(2:nx-1,2:ny-1) * lastw
            domain%ustar(2:nx-1,2:ny-1) = domain%Vm(2:nx-1,1,2:ny-1) * currw
            domain%v10(2:nx-1,2:ny-1) = domain%ustar_new(2:nx-1,2:ny-1) * lastw
            !now calculate master ustar based on U and V combined in quadrature
            domain%wspd3d(2:nx-1,1:nz,2:ny-1) = sqrt(domain%Um(2:nx-1,1:nz,2:ny-1)**2 & 
                                                + domain%Vm(2:nx-1,1:nz,2:ny-1)**2) 
            ! added by Patrik Bohlinger in case we need this later for YSU (some variables seem to 
            ! be 3D in articles like the YSU paper Hong et al. 2006)
            domain%wspd(2:nx-1,2:ny-1) = sqrt(domain%Um(2:nx-1,1,2:ny-1)**2 & 
                                       + domain%Vm(2:nx-1,1,2:ny-1)**2) 
            ! added by Patrik Bohlinger since we need this as input for YSU
             domain%ustar(2:nx-1,2:ny-1) = domain%wspd(2:nx-1,2:ny-1) * currw
            ! ----- end temporary solution ----- !

            ! compute z above ground used for estimating indices for 
            domain%z_agl(2:nx-1,2:ny-1) = (domain%z(2:nx-1,1,2:ny-1)-domain%terrain(2:nx-1,2:ny-1)) 
            !added by Patrik in case we need this later for YSU
            ! calculate the Bulk-Richardson number Rib
            domain%thv(2:nx-1,2:ny-1) = domain%th(2:nx-1,1,2:ny-1) & 
                                        *(1+0.608*domain%qv(2:nx-1,1,2:ny-1)*1000)    
            ! should domain%qv be multiplied by 1000? Did it since domain%qv is in kg/kg and not in g/kg 
            ! normally should be specific humidity and not mixing ratio domain%qv but for first order approach does not matter
            domain%thv3d(2:nx-1,1:nz,2:ny-1) = domain%th(2:nx-1,1:nz,2:ny-1) & 
                                               * (1+0.608*domain%qv(2:nx-1,1:nz,2:ny-1)*1000)   !thv 3D
            domain%thvg(2:nx-1,2:ny-1) = (domain%t2m(2:nx-1,2:ny-1)/domain%pii(2:nx-1,1,2:ny-1)) &
                                         *(1+0.608*domain%qv(2:nx-1,1,2:ny-1)*1000)
            ! t2m should rather be used than skin_t
            domain%thg(2:nx-1,2:ny-1) = domain%t2m(2:nx-1,2:ny-1)/domain%pii(2:nx-1,1,2:ny-1) 
            ! t2m should rather be used than skin_t

            ! variables described for YSU but probably not needed to be calculated outside of the scheme:
            !domain%wstar(2:nx-1,2:ny-1) = domain%ustar(2:nx-1,2:ny-1) / domain%psim(2:nx-1,2:ny-1) ! wstar = vertical wind speed scale
            !domain%thT(2:nx-1,2:ny-1) = propfact * (virtual heat flux)/ domain%wstar ! virtual temperature excess
            !domain%thvg(2:nx-1,2:ny-1) = domain%thv(2:nx-1,2:ny-1) ! for init thvg=thv since thT = 0, 
            !t2m should rather be used than skin_t, b=proportionality factor=7.8, Hong et al, 2006

            ! find value of pbl heights for wspd3d
            !domain%PBLh(2:nx-1,2:ny-1) = Rib_cr * domain%thv(2:nx-1,2:ny-1) * domain%wspd(2:nx-1,2:ny-1)**2 & 
            !                             / gravity * (domain%thv(2:nx-1,2:ny-1) - domain%thvg(2:nx-1,2:ny-1)) !U^2 and thv are from height PBLh in equation

            !write(*,*) "max min domain%pbl_height: ", maxval(domain%pbl_height), minval(domain%pbl_height)
            !write(*,*) "max min domain%PBLh: ", maxval(domain%PBLh), minval(domain%PBLh) ! introduced the 
            !PBLh variabel to not overwrite pbl_height and compare new with old calculations as the pbl 
            !height is one of the most crucial factors of the non-local surface layer calculations needed by the YSU-scheme

            ! Constraint to prevent Rib from becoming too high a lower limit of 0.1 is
            ! applied in for the original surface layer formulation Jiminez et al 2012
            where(domain%wspd(2:nx-1,2:ny-1) < 0.1)
                domain%wspd(2:nx-1,2:ny-1) = 0.1
            endwhere

            domain%Rib(2:nx-1,2:ny-1) = gravity/domain%th(2:nx-1,1,2:ny-1) * domain%z_agl(2:nx-1,2:ny-1) &
                                        * (domain%thv(2:nx-1,2:ny-1) - domain%thvg(2:nx-1,2:ny-1)) &
                                        / domain%wspd(2:nx-1,2:ny-1)**2
            ! From Jiminez et al. 2012, from what height should the theta variables really be, Rib is a function of height? To my understanding the appropriate height is the lower most level.

            ! calculate the integrated similarity functions
            where(domain%Rib(2:nx-1,2:ny-1) >= 0.2)
                !regime = 1, very stable night time conditions
                domain%psim(2:nx-1,2:ny-1) = -10*log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
                domain%psim10(2:nx-1,2:ny-1) = -10*log(10/domain%znt(2:nx-1,2:ny-1))
                domain%psim2m(2:nx-1,2:ny-1) = -10*log(2/domain%znt(2:nx-1,2:ny-1))
                domain%psih(2:nx-1,2:ny-1) = domain%psim(2:nx-1,2:ny-1)
                domain%psih2m(2:nx-1,2:ny-1) = domain%psim2m(2:nx-1,2:ny-1)
                !impose constraints
                where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psim(2:nx-1,2:ny-1) < -10.)
                    domain%psim(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psih(2:nx-1,2:ny-1) < -10.)
                    domain%psih(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psim10(2:nx-1,2:ny-1) < -10.)
                    domain%psim10(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psih2m(2:nx-1,2:ny-1) < -10.)
                    domain%psih2m(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) >= 0.2 .and. domain%psim2m(2:nx-1,2:ny-1) < -10.)
                    domain%psim2m(2:nx-1,2:ny-1) = -10.
                endwhere
            elsewhere (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. domain%Rib(2:nx-1,2:ny-1) > 0.0)
                !regime = 2, damped mechanical turbulence
                domain%psim(2:nx-1,2:ny-1) = -5*domain%Rib(2:nx-1,2:ny-1)*log(domain%z_agl(2:nx-1,2:ny-1) & 
                                             /domain%znt(2:nx-1,2:ny-1))/(1.1-5*domain%Rib(2:nx-1,2:ny-1))
                domain%psim10(2:nx-1,2:ny-1) = -5*domain%Rib(2:nx-1,2:ny-1)*log(10/domain%znt(2:nx-1,2:ny-1)) &
                                             /(1.1-5*domain%Rib(2:nx-1,2:ny-1)) ! Should maybe compute Rib at 10m as well?
                domain%psim2m(2:nx-1,2:ny-1) = -5*domain%Rib(2:nx-1,2:ny-1)*log(2/domain%znt(2:nx-1,2:ny-1)) &
                                             /(1.1-5*domain%Rib(2:nx-1,2:ny-1)) ! Should maybe compute Rib at 2m as well?
                domain%psih(2:nx-1,2:ny-1) = domain%psim(2:nx-1,2:ny-1)
                domain%psih2m(2:nx-1,2:ny-1) = domain%psim2m(2:nx-1,2:ny-1)
                !impose constraints
                where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
                        domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
                        domain%psim(2:nx-1,2:ny-1) < -10.)
                    domain%psim(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
                        domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
                        domain%psih(2:nx-1,2:ny-1) < -10.)
                    domain%psih(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
                        domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
                        domain%psim10(2:nx-1,2:ny-1) < -10.)
                    domain%psim10(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
                        domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
                        domain%psih2m(2:nx-1,2:ny-1) < -10.)
                    domain%psih2m(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0.2 .and. &
                        domain%Rib(2:nx-1,2:ny-1) > 0.0 .and. &
                        domain%psim2m(2:nx-1,2:ny-1) < -10.)
                    domain%psim2m(2:nx-1,2:ny-1) = -10.
                endwhere
            elsewhere (domain%Rib(2:nx-1,2:ny-1).eq.0.0)
                !regime = 3, forced convection
                domain%psim(2:nx-1,2:ny-1) = 0.0
                domain%psim10(2:nx-1,2:ny-1) = 0.0
                domain%psih(2:nx-1,2:ny-1) = 0.0
                domain%psih2m(2:nx-1,2:ny-1) = 0.0
            elsewhere (domain%Rib(2:nx-1,2:ny-1) < 0)
                !regime = 4, free convection
                ! constraints
                where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol(2:nx-1,2:ny-1) < -10.)
                    domain%zol(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol(2:nx-1,2:ny-1) > 0.)
                    domain%zol(2:nx-1,2:ny-1) = 0.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol10(2:nx-1,2:ny-1) < -10.)
                    domain%zol10(2:nx-1,2:ny-1) = -10.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol10(2:nx-1,2:ny-1) > 0.)
                    domain%zol10(2:nx-1,2:ny-1) = 0.
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol2m(2:nx-1,2:ny-1) < -10.)
                    domain%zol2m(2:nx-1,2:ny-1) = -10.              
                endwhere
                where (domain%Rib(2:nx-1,2:ny-1) < 0. .and. domain%zol2m(2:nx-1,2:ny-1) > 0.)
                    domain%zol2m(2:nx-1,2:ny-1) = 0.
                endwhere
                domain%psix(2:nx-1,2:ny-1) = (1.-16.*(domain%zol(2:nx-1,2:ny-1)))**0.25
                domain%psix10(2:nx-1,2:ny-1) = (1.-16.*(domain%zol10(2:nx-1,2:ny-1)))**0.25
                domain%psix2m(2:nx-1,2:ny-1) = (1.-16.*(domain%zol2m(2:nx-1,2:ny-1)))**0.25
                domain%psim(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix(2:nx-1,2:ny-1))/2.) & 
                                             + log((1.+domain%psix(2:nx-1,2:ny-1)**2.)/2.) &
                                             - 2.*atan(domain%psix(2:nx-1,2:ny-1))+pi/2.
                domain%psim10(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix10(2:nx-1,2:ny-1))/2.) &
                                             + log((1.+domain%psix10(2:nx-1,2:ny-1)**2.)/2.) &
                                             - 2.*atan(domain%psix10(2:nx-1,2:ny-1))+pi/2.
                domain%psih(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix(2:nx-1,2:ny-1)**2.)/2.)
                domain%psih2m(2:nx-1,2:ny-1) = 2.*log((1.+domain%psix2m(2:nx-1,2:ny-1)**2.)/2.)
            endwhere
            ! constrain psim and psih also for 2m and 10m
            where (domain%psim(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
                domain%psim(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
            endwhere
            where (domain%psim2m(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
                domain%psim2m(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
            endwhere
            where (domain%psim10(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
                domain%psim10(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
            endwhere
            where (domain%psih(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
                domain%psih(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
            endwhere
            where (domain%psih2m(2:nx-1,2:ny-1) > 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)))
                domain%psih2m(2:nx-1,2:ny-1) = 0.9 * log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
            endwhere

            ! calculate thstar = temperature scale
            domain%thstar(2:nx-1,2:ny-1) = karman*(domain%th(2:nx-1,1,2:ny-1)-domain%thg) &
                                           / log(domain%z_agl(2:nx-1,2:ny-1) & 
                                           / domain%znt(2:nx-1,2:ny-1))-domain%psih(2:nx-1,2:ny-1)
            ! Averaging ustar with ustar from previous time step to suppress
            ! large oscillations
            domain%ustar_tmp(2:nx-1,2:ny-1) = karman*domain%wspd(2:nx-1,2:ny-1)/(log(domain%z_agl(2:nx-1,2:ny-1) & 
                                              / domain%znt(2:nx-1,2:ny-1))-domain%psim(2:nx-1,2:ny-1))
            domain%ustar_new(2:nx-1,2:ny-1) = (domain%ustar_tmp(2:nx-1,2:ny-1) + domain%ustar_new(2:nx-1,2:ny-1))/2.0
            ! preventing ustar from being smaller than 0.1 as it could be under
            ! very stable conditions, Jiminez et al. 2012
            where(domain%ustar_new(2:nx-1,2:ny-1) < 0.1)
                domain%ustar_new(2:nx-1,2:ny-1) = 0.1
            endwhere
            ! calculate the Monin-Obukhov  stability parameter zol (z over l)
            ! using ustar from the similarity theory
            domain%zol(2:nx-1,2:ny-1) = (karman*gravity*domain%z_agl(2:nx-1,2:ny-1))/domain%th(2:nx-1,1,2:ny-1) &
                                        * domain%thstar(2:nx-1,2:ny-1)/(domain%ustar_new(2:nx-1,2:ny-1) &
                                        * domain%ustar_new(2:nx-1,2:ny-1))
            domain%zol10(2:nx-1,2:ny-1) = (karman*gravity*10)/domain%th(2:nx-1,1,2:ny-1) * domain%thstar(2:nx-1,2:ny-1) &
                                        / (domain%ustar_new(2:nx-1,2:ny-1)*domain%ustar_new(2:nx-1,2:ny-1))
            domain%zol2m(2:nx-1,2:ny-1) = (karman*gravity*2)/domain%th(2:nx-1,1,2:ny-1) * domain%thstar(2:nx-1,2:ny-1) &
                                        / (domain%ustar_new(2:nx-1,2:ny-1)*domain%ustar_new(2:nx-1,2:ny-1))
            ! calculate pblh over l using ustar and thstar from the similarity theory
            domain%hol(2:nx-1,2:ny-1) = (karman*gravity*domain%PBLh(2:nx-1,2:ny-1))/domain%th(2:nx-1,1,2:ny-1) &
                                        * domain%thstar(2:nx-1,2:ny-1)/(domain%ustar_new(2:nx-1,2:ny-1) &
                                        * domain%ustar_new(2:nx-1,2:ny-1))
            ! arbitrary variables
            domain%gz1oz0(2:nx-1,2:ny-1)=log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))
            ! calculating dtmin
            domain%dtmin = domain%dt / 60.0
            ! p_top as a scalar, choosing just minimum from ptop as a start
            p_top = minval(domain%ptop)
            ! compute the dimensionless bulk coefficent for momentum, heat and
            ! moisture
            !domain%exch_m(2:nx-1,2:ny-1) = (karman**2)/(log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)) &
            !- domain%psim(2:nx-1,2:ny-1))**2
            !domain%exch_h(2:nx-1,2:ny-1) = (karman**2)/((log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))&
            !- domain%psim(2:nx-1,2:ny-1))*(log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1))-domain%psih(2:nx-1,2:ny-1)))
            !domain%exch_q(2:nx-1,2:ny-1) = (karman**2)/((log(domain%z_agl(2:nx-1,2:ny-1)/domain%znt(2:nx-1,2:ny-1)) &
            !- domain%psim(2:nx-1,2:ny-1)) * (log(domain%rho(2:nx-1,1,2:ny-1)*cp*karman*domain%ustar_new(2:nx-1,2:ny-1) &
            !*domain%z_agl(2:nx-1,2:ny-1)/cs)-psih(2:nx-1,2:ny-1)))

            ! counter is just a variable helping me to detect how much rounds this subroutine went through
            write(*,*) "Counter: ", counter
            counter = counter + 1
            write(*,*) "Counter: ", counter
            ! end surface layer calculations introduced by Patrik Bohlinger
        endif        
        write(*,*) "end surface layer calculations"
        ! ----- end sfc layer calculations ----- !

        ! finally, calculate the real vertical motions (including U*dzdx + V*dzdy)
        lastw=0
        do z=1,nz
            ! compute the U * dz/dx component of vertical motion
            uw    = domain%u(2:nx,  z,2:ny-1) * domain%dzdx(:,2:ny-1)
            ! compute the V * dz/dy component of vertical motion
            vw    = domain%v(2:nx-1,z,2:ny  ) * domain%dzdy(2:nx-1,:)
            ! convert the W grid relative motion to m/s
            currw = domain%w(2:nx-1,z,2:ny-1) * domain%dz_inter(2:nx-1,z,2:ny-1) / domain%dx
            
            ! compute the real vertical velocity of air by combining the different components onto the mass grid
            domain%w_real(2:nx-1,z,2:ny-1) = (uw(1:nx-2,:) + uw(2:nx-1,:))*0.5 &
                                            +(vw(:,1:ny-2) + vw(:,2:ny-1))*0.5 &
                                            +(lastw + currw) * 0.5
                                            
            lastw=currw ! could avoid this memcopy cost using pointers or a single manual loop unroll
        end do
        
    end subroutine diagnostic_update


    !>------------------------------------------------------------
    !!  Uses the dt from step() to convert the forcing increments to per/time step increments
    !!
    !!  Divides dXdt variables by n timesteps so that after adding it N times we will be at the
    !!  correct final value. 
    !!
    !! @param bc        Boundary conditions structure containing dXdt variables
    !! @param nsteps    Number of time steps calculated in step() function
    !! @param options   Model options structure specifies which variables need to be updated
    !!
    !!------------------------------------------------------------
    subroutine apply_dt(bc,nsteps, options)
        implicit none
        type(bc_type), intent(inout) :: bc
        integer,intent(in)::nsteps
        type(options_type), intent(in) :: options
        integer::j, ny, nx
        ny=size(bc%du_dt,3)
        nx=size(bc%du_dt,1)

        bc%du_dt  = bc%du_dt  / nsteps
        bc%dv_dt  = bc%dv_dt  / nsteps
        bc%dp_dt  = bc%dp_dt  / nsteps
        bc%dth_dt = bc%dth_dt / nsteps
        bc%dqv_dt = bc%dqv_dt / nsteps
        bc%dqc_dt = bc%dqc_dt / nsteps

        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%landsurface==kLSM_BASIC) then
            bc%dsh_dt   = bc%dsh_dt   / nsteps
            bc%dlh_dt   = bc%dlh_dt   / nsteps
        endif
        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%boundarylayer==kPBL_BASIC) then
            bc%dpblh_dt = bc%dpblh_dt / nsteps
        endif
        ! these only get updated if we are using the fluxes derived from the forcing model
        if (options%physics%radiation==kRA_BASIC) then
            bc%dsw_dt   = bc%dsw_dt   / nsteps
            bc%dlw_dt   = bc%dlw_dt   / nsteps
        endif
        bc%dsst_dt   = bc%dsst_dt   / nsteps

        ! parallel version didn't seem to work
!         !$omp parallel firstprivate(ny,nsteps) &
!         !$omp private(j) &
!         !$omp shared(bc)
!         !$omp do schedule(static)
!         do j=1,ny
!             bc%du_dt(:,:,j)  = bc%du_dt(:,:,j)  / nsteps
!             bc%dv_dt(:,:,j)  = bc%dv_dt(:,:,j)  / nsteps
!             bc%dp_dt(:,:,j)  = bc%dp_dt(:,:,j)  / nsteps
!             bc%dth_dt(:,j,:) = bc%dth_dt(:,j,:) / nsteps
!             bc%dqv_dt(:,j,:) = bc%dqv_dt(:,j,:) / nsteps
!             bc%dqc_dt(:,j,:) = bc%dqc_dt(:,j,:) / nsteps
!
!             ! these only get updated if we are using the fluxes derived from the forcing model
!             if (options%physics%landsurface==kLSM_BASIC) then
!                 bc%dsh_dt(:,j)   = bc%dsh_dt(:,j)   / nsteps
!                 bc%dlh_dt(:,j)   = bc%dlh_dt(:,j)   / nsteps
!             endif
!             ! these only get updated if we are using the fluxes derived from the forcing model
!             if (options%physics%boundarylayer==kPBL_BASIC) then
!                 bc%dpblh_dt(:,j) = bc%dpblh_dt(:,j) / nsteps
!             endif
!             ! these only get updated if we are using the fluxes derived from the forcing model
!             if (options%physics%radiation==kRA_BASIC) then
!                 bc%dsw_dt(:,j)   = bc%dsw_dt(:,j)   / nsteps
!                 bc%dlw_dt(:,j)   = bc%dlw_dt(:,j)   / nsteps
!             endif
!             bc%dsw_dt(:,j)   = bc%dsw_dt(:,j)       / nsteps
!             
!         end do
!         !$omp end do
!         !$omp end parallel
!         if (nx>ny) then
!             do j=ny+1,nx
!                 bc%dth_dt(:,j,:) = bc%dth_dt(:,j,:) / nsteps
!                 bc%dqv_dt(:,j,:) = bc%dqv_dt(:,j,:) / nsteps
!                 bc%dqc_dt(:,j,:) = bc%dqc_dt(:,j,:) / nsteps
!             end do
!         endif
            
    end subroutine apply_dt

    !>------------------------------------------------------------
    !!  Calculate the maximum stable time step given some CFL criteria 
    !!  
    !!  For each grid cell, find the mean of the wind speeds from each
    !!  direction * sqrt(3) for the 3D advection CFL limited time step
    !!  Also find the maximum wind speed anywhere in the domain to check 
    !!  against a 1D advection limit. 
    !!  
    !! @param dx  [ scalar ]        horizontal grid cell width  [m]
    !! @param u   [nx+1 x nz x ny]  east west wind speeds       [m/s]
    !! @param v   [nx x nz x ny+1]  North South wind speed      [m/s]
    !! @param w   [nx x nz x ny]    vertical wind speed         [m/s]
    !! @param CFL [ scalar ]        CFL limit to use (e.g. 1.0)
    !! @return dt [ scalar ]        Maximum stable time step    [s]
    !!
    !!------------------------------------------------------------
    function compute_dt(dx, u, v, w, CFL, cfl_strictness) result(dt)
        real, intent(in) :: dx
        real, intent(in), dimension(:,:,:) :: u,v,w
        real, intent(in) :: CFL
        integer, intent(in) :: cfl_strictness
        ! output value
        real :: dt
        ! locals
        real :: three_d_cfl = 0.577350269 ! = sqrt(3)/3
        integer :: i, j, k, nx, nz, ny, zoffset
        real :: maxwind3d, maxwind1d, current_wind, sqrt3
        
        sqrt3 = sqrt(3.0) * 1.001 ! with a safety factor
        
        maxwind1d = 0
        maxwind3d = 0
        
        nx = size(w,1)
        nz = size(w,2)
        ny = size(w,3)

        if (cfl_strictness==1) then
            ! to ensure we are stable for 1D advection:
            maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
            maxwind1d = max( maxwind1d, maxval(abs(w)))
            
            maxwind3d = maxwind1d * sqrt3
        else if (cfl_strictness==5) then
            maxwind3d = maxval(abs(u)) + maxval(abs(v)) + maxval(abs(w))
        else
            ! to ensure we are stable for 3D advection we'll use the average "max" wind speed
            ! but that average has to be divided by sqrt(3) for stability in 3 dimensional advection
            do j=1,ny
                do k=1,nz
                    if (k==1) then
                        zoffset = 0
                    else
                        zoffset = -1
                    endif
                    
                    do i=1,nx
                        ! just compute the sum of the wind speeds, but take the max over the two
                        ! faces of the grid cell (e.g. east and west sides)
                        ! this will be divided by 3 later by three_d_cfl
                        current_wind = max(abs(u(i,k,j)), abs(u(i+1,k,j))) &
                                      +max(abs(v(i,k,j)), abs(v(i,k,j+1))) & 
                                      +max(abs(w(i,k,j)), abs(w(i,k+zoffset,j)))
                        maxwind3d = max(maxwind3d, current_wind)
                    ENDDO
                ENDDO
            ENDDO
            
            if (cfl_strictness==2) then
                ! effectively divides by 3 to take the mean and multiplies by the sqrt(3) for the 3D advection limit
                maxwind3d = maxwind3d * three_d_cfl
                
                ! to ensure we are stable for 1D advection:
                maxwind1d = max( maxval(abs(u)), maxval(abs(v)))
                maxwind1d = max( maxwind1d, maxval(abs(w)))
                ! also insure stability for 1D advection
                maxwind3d = max(maxwind1d,maxwind3d)

            ! else if (cfl_strictness==3) then 
            !   leave maxwind3d as the sum of the max winds
            ! This should be the default, does it need to be multiplied by sqrt(3)?
            elseif (cfl_strictness==4) then
                maxwind3d = maxwind3d * sqrt3
            endif
                
        endif        
        
        dt = CFL * dx / maxwind3d
        
    end function compute_dt
    
    !>------------------------------------------------------------
    !!  Step forward one IO time step. 
    !! 
    !!  Calculated the internal model time step to satisfy the CFL criteria, 
    !!  then updates all forcing update increments for that dt and loops through
    !!  time calling physics modules. 
    !!  Also checks to see if it is time to write a model output file.
    !!
    !! @param domain    domain data structure containing model state
    !! @param options   model options structure
    !! @param bc        model boundary conditions data structure
    !! @param model_time    Current internal model time (seconds since start of run)
    !! @param next_output   Next time to write an output file (in "model_time")
    !!
    !!------------------------------------------------------------
    subroutine step(domain,options,bc,model_time,next_output)
        implicit none
        type(domain_type),intent(inout)::domain
        type(bc_type),intent(inout)::bc
        type(options_type),intent(in)::options
        real*8,intent(inout)::model_time,next_output
        real*8::end_time
        integer::i,ntimesteps,tenp
        real::dt,dtnext, CFL, cfl_reduction_factor
        
        CFL = 1.0
        cfl_reduction_factor = options%cfl_reduction_factor
        CFL = CFL * cfl_reduction_factor
        
!       compute internal timestep dt to maintain stability
!       courant condition for 3D advection. Note that w is normalized by dx/dz
!       pick the minimum dt from the begining or the end of the current timestep
        dt = compute_dt(domain%dx, domain%u, domain%v, domain%w, CFL, cfl_strictness=options%cfl_strictness)
        ! print*, "Current dt = ",dt
        dtnext = compute_dt(domain%dx, bc%next_domain%u, bc%next_domain%v, bc%next_domain%w, CFL, cfl_strictness=options%cfl_strictness)
        ! print*, "Next    dt = ",dtnext
        dt=min(dt,dtnext)
!       set an upper bound on dt to keep microphysics and convection stable (?) not sure what time is required here. 
        dt=min(dt,180.0) !better min=180?
!       if we have too small a time step just throw an error
        if (dt<1e-1) then
            write(*,*) "dt=",dt
            write(*,*) "Umax",maxval(abs(domain%u)),maxval(abs(bc%next_domain%u))
            write(*,*) "Vmax",maxval(abs(domain%v)),maxval(abs(bc%next_domain%v))
            write(*,*) "Wmax",maxval(abs(domain%w)),maxval(abs(bc%next_domain%w))
            call write_domain(domain,options,99998)
            call write_domain(bc%next_domain,options,99999,'timestep_error_file.nc')
            stop "ERROR time step too small"
        endif
        
!       make dt an integer fraction of the full timestep
        dt = options%in_dt/ceiling(options%in_dt/dt)
!       calculate the number of timesteps
        ntimesteps=nint(options%in_dt/dt)
        end_time=model_time+options%in_dt
        
!       adjust the boundary condition dXdt values for the number of time steps
        call apply_dt(bc,ntimesteps,options)
        write(*,*) "    dt=",dt, "nsteps=",ntimesteps
        
!       ensure internal model consistency (should only need to be called here when the model starts...)
!       e.g. for potential temperature and temperature
        call diagnostic_update(domain,options)
        
!       now just loop over internal timesteps computing all physics in order (operator splitting...)
        do i=1,ntimesteps
            model_time=model_time+dt
            if (dt>1e-5) then
                write(*,*) "--- before advect ---"
                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib),MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol),MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol),MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt),MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z) 
                write(*,*) "domain%z_agl: ",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th),MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv),MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl),MINVAL(domain%tend%qv_pbl)
                write(*,*) "--- call advect ---"
                call advect(domain,options,dt)
                write(*,*) "--- after advect ---"
                write(*,*) "--- before mp ---"
                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
                write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th),MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv),MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl),MINVAL(domain%tend%qv_pbl)
                write(*,*) "--- call mp ---"
                call mp(domain,options,dt, model_time)
                write(*,*) "--- after mp ---"
                write(*,*) "--- before rad ---"
                call rad(domain,options,model_time/86400.0+50000, dt)
                write(*,*) "--- before lsm ---"
                call lsm(domain,options,dt,model_time)
                write(*,*) "--- before pbl ---"
                !call io_write("bc%dqv_dt_bpbl.nc","data",bc%dqv_dt)
                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
                write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th),MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv),MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl),MINVAL(domain%tend%qv_pbl)
                write(*,*) "--- call pbl ---"
                call pbl(domain,options,dt)
                write(*,*) "--- after pbl ---"
                !call io_write("bc%dqv_dt_apbl.nc","data",bc%dqv_dt)
                write(*,*) "--- before convect ---"
                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
                write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th),MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv),MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl),MINVAL(domain%tend%qv_pbl)
                write(*,*) "--- call convect ---"
                call convect(domain,options,dt)
                write(*,*) "--- after convect ---"
                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
                write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th),MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv),MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl),MINVAL(domain%tend%qv_pbl)
    !           apply/update boundary conditions including internal wind and pressure changes. 
                write(*,*) "--- before forcing_update ---"
                !call io_write("domain%exch_h_bfu.nc","data",domain%exch_h)
                !call io_write("domain%th_bfu.nc","data",domain%th)
                !call io_write("domain%cloud_bfu.nc","data",domain%cloud)
                !call io_write("domain%qv_bfu.nc","data",domain%qv)
                call io_write("domain%PBLh_bfu.nc","data",domain%PBLh)
                !call io_write("bc%dqv_dt_bfu.nc","data",bc%dqv_dt)

                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
                write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th),MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv),MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl),MINVAL(domain%tend%qv_pbl)
                write(*,*) "--- call forcing_update ---"
                call forcing_update(domain,bc,options)
                write(*,*) "--- after forcing update ---"
                !call io_write("domain%exch_h_afu.nc","data",domain%exch_h)
                !call io_write("domain%th_afu.nc","data",domain%th)
                !call io_write("domain%cloud_afu.nc","data",domain%cloud)
                !call io_write("domain%qv_afu.nc","data",domain%qv)
                call io_write("domain%PBLh_afu.nc","data",domain%PBLh)

                write(*,*) "domain%t: ", MAXVAL(domain%t), MINVAL(domain%t)
                write(*,*) "domain%th: ", MAXVAL(domain%th), MINVAL(domain%th)
                write(*,*) "domain%Um: ", MAXVAL(domain%Um), MINVAL(domain%Um)
                write(*,*) "domain%Vm: ", MAXVAL(domain%Vm), MINVAL(domain%Vm)
                write(*,*) "domain%u: ", MAXVAL(domain%u), MINVAL(domain%u)
                write(*,*) "domain%v: ", MAXVAL(domain%v), MINVAL(domain%v)
                write(*,*) "domain%qv: ", MAXVAL(domain%qv), MINVAL(domain%qv)
                write(*,*) "domain%cloud: ", MAXVAL(domain%cloud),MINVAL(domain%cloud)
                write(*,*) "domain%ice: ", MAXVAL(domain%ice), MINVAL(domain%ice)
                write(*,*) "domain%psim: ", MAXVAL(domain%psim),MINVAL(domain%psim)
                write(*,*) "domain%psih: ", MAXVAL(domain%psih),MINVAL(domain%psih)
                write(*,*) "domain%PBLh: ", MAXVAL(domain%PBLh),MINVAL(domain%PBLh)
                write(*,*) "domain%Rib: ", MAXVAL(domain%Rib), MINVAL(domain%Rib)
                write(*,*) "domain%hol: ", MAXVAL(domain%hol), MINVAL(domain%hol)
                write(*,*) "domain%zol: ", MAXVAL(domain%zol), MINVAL(domain%zol)
                write(*,*) "domain%znt: ", MAXVAL(domain%znt), MINVAL(domain%znt)
                write(*,*) "domain%ustar: ",MAXVAL(domain%ustar),MINVAL(domain%ustar)
                write(*,*) "domain%ustar_new: ",MAXVAL(domain%ustar_new),MINVAL(domain%ustar_new)
                write(*,*) "domain%exch_h: ",MAXVAL(domain%exch_h),MINVAL(domain%exch_h)
                write(*,*) "domain%z: ", MAXVAL(domain%z),MINVAL(domain%z)
                write(*,*) "domain%z_agl:",MAXVAL(domain%z_agl),MINVAL(domain%z_agl)
                write(*,*) "domain%tend%th: ", MAXVAL(domain%tend%th), MINVAL(domain%tend%th)
                write(*,*) "domain%tend%qv: ", MAXVAL(domain%tend%qv), MINVAL(domain%tend%qv)
                write(*,*) "domain%tend%qv_pbl: ", MAXVAL(domain%tend%qv_pbl), MINVAL(domain%tend%qv_pbl)
                write(*,*) "done, I guess ..."
                !call io_write("pbl_qv_tendency.nc","data",domain%tend%qv_pbl)
            
    !           step model time forward
                domain%model_time=model_time
                if ((abs(model_time-next_output)<1e-1).or.(model_time>next_output)) then
                    call write_domain(domain,options,nint((model_time-options%time_zero)/options%out_dt))
                    next_output=next_output+options%out_dt
                endif
    !           in case out_dt and in_dt arent even multiples of each other.  Make sure we don't over step
                if ((model_time+dt)>end_time) then
                    dt=end_time-model_time
                endif
            endif
        enddo
        
    end subroutine step
end module time_step
