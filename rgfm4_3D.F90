subroutine RGFM(mype,nx,ny,nz,curvature,steps)    

  !    subroutine RGFM(mype,n1,ny,vnorm,steps,irk)

  use paramesh_dimensions
  use physicaldata
  use tree
  use workspace
  use paramesh_interfaces

  include 'mpif.h'

  integer :: mype,i,j,k,ii,jj,nrpocs,ierr,mark(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),nLayer,&
       markStripP(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),markStripM(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       iopt,nlayers,lb,l,m,n
  integer,intent(in) :: steps

  Real :: p,q,r,dx,dy,dz,w1,w2,lS

  real :: dens_oA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),press_oA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       dens_oB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),press_oB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       vx_oA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vy_oA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       vx_oB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vy_oB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       vz_oA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vz_oB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       vn_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),   vn_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       vof(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)
  real :: dens_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vx_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       dens_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vx_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &                   
       press_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vy_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       press_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vy_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       dens_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vx_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       press_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vy_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vz_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vz_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),vz_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)


  real :: dens_IL,vn_I,press_IL,dens_IR,press_IR,&
       xI,yI,zI,dens_A,vx_A,vy_A,vz_A,vn_A,press_A, &
       dens_B,vx_B,vy_B,vz_B,vn_B,press_B,E_cell,hStrip
  real,intent(in) ::  nx(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       ny(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),nz(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       curvature(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)


  real :: vn_IL,vn_IR,kappa,ph_,ph_ip1,ph_jp1,ph_kp1    
  real :: rho,rhoVelX1,rhoVelY1,rho1E1


  hStrip=1.0 

  kappa=0. 
  nLayer=0

  !((((DEBUGGING
  if (mype==0) print *, 'entering RGFM code'
  !))))) END DEBUGGING

!!!!****note GA, GB subscripts do not mean values actually in ghost regions, until these GA, GB vars. GA is originally the
  !values of U_A and 6dx into U_B. They will have the last node before the interface modified by Riemann solver.
  !These are extrapolated across the interface much later in the code.

  !    print *, '------------'

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  mark = 0
  markStripP=0
  markStripM=0

  dens_GA = 0.  ; dens_GB = 0.
  vx_GA = 0.0    ; vx_GB = 0.0
  vy_GA = 0.0    ; vy_GB = 0.0
  vz_GA = 0.0    ; vz_GB = 0.0
  press_GA = 0.0 ; press_GB = 0.0

  dens_new=0.0;
  press_new=0.0;
  vx_new=0.0;
  vy_new=0.0;
  vz_new=0.0;
  vof=0.0;
  !vofP=0.0;vofM=0.0;
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  iopt=1
  nguard = 3
  nlayers=nguard

  call amr_guardcell(mype,iopt,nlayers)

  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then
        do i=1,nxb+nguard+2
           do j=1,nyb+nguard+2
              do k=1,nzb+nguard+2
                 if ( abs(unk(11,i,j,k,lb)) < 2.598*dx_min ) then
!!$                 do l=-1,1
!!$                    do m=-1,1
!!$                       do n = -1,1
!!$                          
!!$                          if ( unk(11,i,j,k,lb)*unk(11,i+l,j+m,k+n,lb) .le. 0.0 ) then !this edits only cells bordering sign change in phi, including diagonals
!!$                             mark(i,j,k,lb) = 1
!!$                          endif
!!$                       enddo
!!$                    enddo
!!$                 enddo


                    if (unk(11,i+1,j,k,lb)*unk(11,i,j,k,lb) <=0. ) then
                       mark(i,j,k,lb) = 1
                       mark(i+1,j,k,lb) = 1
                    endif
                    if (unk(11,i,j+1,k,lb)*unk(11,i,j,k,lb) <=0. ) then
                       mark(i,j,k,lb) = 1
                       mark(i,j+1,k,lb) = 1
                    endif

                    if (unk(11,i,j,k+1,lb)*unk(11,i,j,k,lb) <=0. ) then
                       mark(i,j,k,lb) = 1
                       mark(i,j,k+1,lb) = 1
                    endif



                 endif


              enddo
           enddo
        enddo
     endif
  enddo



  if ( mype == 0 ) then 
     print *, 'exiting vof calc'
  endif
  call vof_calc(nx,ny,nz,vof) !calculating volume fraction 1dx away from 0 level set in - dir
  !note the unk are already smoothed in the IC
  t_unk=unk !store temp copy so we can use paramesh machinery on unk
  !    lSeg=0.0


  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then
        do k=1,nzb+2*nguard
           do j=1,nyb+2*nguard              
              do i=1,nxb+2*nguard
                 !if ( abs(unk(11,i,j,k,lb)) <= 4.0*dx_min ) then
                 !if (mark(i,j,k,lb) = 1)
                 if (unk(11,i,j,k,lb)<=4.0*dx_min) then

                    !dens_GA((i-nguard):(i+nguard),(j-nguard):j+nguard),(k-nguard):(k+nguard),lb)  = unk(1,(i-nguard):(i+nguard),(j-nguard):j+nguard),(k-nguard):(k+nguard),lb)
                    dens_GA(i,j,k,lb) = unk(1,i,j,k,lb)
                    vx_GA(i,j,k,lb)    = unk(2,i,j,k,lb)/unk(1,i,j,k,lb)
                    vy_GA(i,j,k,lb)    = unk(3,i,j,k,lb)/unk(1,i,j,k,lb)
                    vz_GA(i,j,k,lb) =  unk(4,i,j,k,lb)/unk(1,i,j,k,lb)


                    press_GA(i,j,k,lb) = (gamm_A-1)*(unk(5,i,j,k,lb)-0.5*(unk(2,i,j,k,lb)**2 + &
                         unk(3,i,j,k,lb)**2+unk(4,i,j,k,lb)**2)/unk(1,i,j,k,lb) )-gamm_A*p_inf_A
                    !endif
                    if (unk(11,i,j,k,lb) < 0.) then
                       if ( press_GA(i,j,k,lb) <0. )then
                          print *, 'press_GA<0 rgfm ln 167'
                          print *, '(i,j,k,lb) = ', i,j,k,lb
                          print *, 'dens_GA(i,j,k,lb) = ',dens_GA(i,j,k,lb)
                          print *, 'vx_GA(i,j,k,lb) = ',vx_GA(i,j,k,lb)
                          print *, 'vy_GA(i,j,k,lb) = ',vy_GA(i,j,k,lb)
                          print *, 'vz_GA(i,j,k,lb) = ', vz_GA(i,j,k,lb)
                          print *, 'gamm_A = ', gamm_A, 'p_inf_A = ', p_inf_A
                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)

                       endif

                       if ( (abs(press_GA(i,j,k,lb) ) > 45098. ) .and.  (abs(press_GA(i,j,k,lb) ) < 45099. ) ) then
                          print *, ' press_GA(i,j,k,lb) = ', press_GA(i,j,k,lb) 
                          print *, '(i,j,k,lb) = ', i,j,k,lb
                          print *, 'dens_GA(i,j,k,lb) = ',dens_GA(i,j,k,lb)
                          print *, 'vx_GA(i,j,k,lb) = ',vx_GA(i,j,k,lb)
                          print *, 'vy_GA(i,j,k,lb) = ',vy_GA(i,j,k,lb)
                          print *, 'vz_GA(i,j,k,lb) = ', vz_GA(i,j,k,lb)
                          print *, 'gamm_A = ', gamm_A, 'p_inf_A = ', p_inf_A

                       endif
                    endif
                    press_oA(i,j,k,lb) = press_GA(i,j,k,lb)
                    dens_oA(i,j,k,lb)=unk(1,i,j,k,lb)
                 endif
                 if (unk(11,i,j,k,lb) >=-4.0*dx_min) then
                    dens_GB(i,j,k,lb)  = unk(6,i,j,k,lb)
                    vx_GB(i,j,k,lb)    = unk(7,i,j,k,lb)/dens_GB(i,j,k,lb)
                    vy_GB(i,j,k,lb)    = unk(8,i,j,k,lb)/dens_GB(i,j,k,lb)
                    vz_GB(i,j,k,lb) = unk(9,i,j,k,lb)/dens_GB(i,j,k,lb)
                    ! vn_GB(i,j,k,lb)    = vx_GB(i,j,k,lb)*nx(i,j,k,lb)+vy_GB(i,j,k,lb)*ny(i,j,k,lb)+vz_GB(i,j,k,lb)*nz(i,j,k,lb) !v_GB dot n_hat
                    !normals are not filled in guardcells!
!!$                       if ( (mark(i,j,k,lb) == 1) .and. (unk(11,i,j,k,lb) > 0.) ) then
!!$                       press_GB(i,j,k,lb) = ((gamm_B*(1.-vof(i,j,k,lb) )+gamm_A*vof(i,j,k,lb)-1.))*&
!!$                            (unk(10,i,j,k,lb)-0.5*(unk(7,i,j,k,lb)**2 + &
!!$                            unk(8,i,j,k,lb)**2+unk(9,i,j,k,lb)**2)/dens_GB(i,j,k,lb)) &
!!$                            -( (1. - vof(i,j,k,lb))*gamm_B*p_inf_B + vof(i,j,k,lb)*gamm_A*p_inf_A )

                    ! else
                    press_GB(i,j,k,lb) = (gamm_B-1)*(unk(10,i,j,k,lb)-0.5*(unk(7,i,j,k,lb)**2 + &
                         unk(8,i,j,k,lb)**2+unk(9,i,j,k,lb)**2)/dens_GB(i,j,k,lb))-gamm_B*p_inf_B
                    !endif
                    if (unk(11,i,j,k,lb) > 0.) then
                       if (press_GB(i,j,k,lb) < 0.01) then
                          print *, 'press_GB(i,j,k,lb) < 0.01, press_GB = ', press_GB(i,j,k,lb), 'RGFM ln 209'
                          print *, 'i,j,k,lb = ', i,j,k,lb
                          
                          print *, 'dens_GB(i,j,k,lb) = ',dens_GB(i,j,k,lb)
                          print *, 'vx_GB(i,j,k,lb) = ',vx_GB(i,j,k,lb)
                          print *, 'vy_GB(i,j,k,lb) = ',vy_GB(i,j,k,lb)
                          print *, 'vz_GB(i,j,k,lb) = ', vz_GB(i,j,k,lb)
                          print *, 'gamm_B = ', gamm_B, 'p_inf_B = ', p_inf_B
                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
                       endif
                    endif
                    press_oB(i,j,k,lb) = press_GB(i,j,k,lb)
                    dens_oB(i,j,k,lb)=unk(6,i,j,k,lb)
                 endif




                 !end if
              end do !k
           end do !j
        end do !i
     end if !nodetype.eq.1
  enddo !lb

  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then
        do i=1+nguard,nxb+nguard
           do j=1+nguard,nyb+nguard
              do k=1+nguard,nzb+nguard

                 if ( mark(i,j,k,lb) ==1) then
                    if (abs(unk(11,i,j,k,lb))>1.2*sqrt(3.)*dx_min ) then
                       mark(i,j,k,lb) = 0
                    else
                       
                    dx = bsize(1,lb)/real(nxb)
                    dy = bsize(2,lb)/real(nyb)
                    dz = bsize(3,lb)/real(nzb)
                    xi = bnd_box(1,1,lb) + dx*(real(i-nguard)-0.5)
                    yi = bnd_box(1,2,lb) + dy*(real(j-nguard)-0.5)        
                    zi = bnd_box(1,3,lb) + dz*(real(k-nguard)-0.5)
                    if (curvature(i,j,k,lb) .eq. 12345678. ) then
                       print *, 'curvature out of bounds RGFM ln 223'
                       print *, '(i,j,k) = (' ,i,j,k, ')'
                       Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
                    endif




                    if(unk(11,i,j,k,lb)<=0.) then

!!$                       ! BEGIN DEBUG
!!$                       if ( press_GA(i,j,k,lb) < 100000. ) then
!!$                          print *, 'press too low rgfm 154' , press_GA(i,j,k,lb)
!!$                       endif
                       if ( dens_GA(i,j,k,lb) - .0001  < .000 ) then
                          print *, 'dens too low rgfm 157' , dens_GA(i,j,k,lb)
                       endif
!!$                       if ( abs(vx_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vx too low rgfm 160' , vx_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vy_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vy too low rgfm 160' , vy_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vz_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vz too low rgfm 160' , vz_GA(i,j,k,lb)
!!$                       endif
!!$                       ! END DEBUG


                       !call steps 1, 2,3 on paper
                       call setUpDataForRiemannProblem(xi,yi,zi,'A',i,j,k,lb,nx,ny,nz,curvature,&
                            dens_GA,press_GA,vx_GA,vy_GA,vz_GA,dens_GB,press_GB,vx_GB,vy_GB,vz_GB,&
                            dens_new,vx_new,vy_new,vz_new,press_new)
!!$                       ! BEGIN DEBUG
!!$                       if ( press_GA(i,j,k,lb) < 100000. ) then
!!$                          print *, 'pressGA too low rgfm 154' , press_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( dens_GA(i,j,k,lb) - 1.  < .000 ) then
!!$                          print *, 'densGA too low rgfm 157' , dens_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vx_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vxGA too low rgfm 160' , vx_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vy_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vyGA too low rgfm 160' , vy_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vz_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vzGA too low rgfm 160' , vz_GA(i,j,k,lb)
!!$                       endif
!!$                       ! END DEBUG
                       ! print *, dens_new(i,j,k,lb),vx_new(i,j,k,lb),vy_new(i,j,k,lb)&
                       !      ,vz_new(i,j,k,lb),press_new(i,j,k,lb),'exiting setUpData in RGFM flA'
                    elseif ( unk(11,i,j,k,lb) >= 0. ) then

!!$                       ! BEGIN DEBUG
!!$                       if ( press_GB(i,j,k,lb) > 10000.1 ) then
!!$                          print *, 'pressGB too low rgfm 181 post riemann' , press_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( dens_GB(i,j,k,lb) - 0.125  < .000 ) then
!!$                          print *, 'densGB too low rgfm 184 post riemann' , dens_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vx_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vxGB too low rgfm 187 post riemann' , vx_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vy_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vyGB too low rgfm 190 post riemann' , vy_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vz_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vzGB too low rgfm 193 post riemann' , vz_GB(i,j,k,lb)
!!$                       endif
!!$                       ! END DEBUG
                       call setUpDataForRiemannProblem(xi,yi,zi,'B',i,j,k,lb,nx,ny,nz,curvature,&
                            dens_GA,press_GA,vx_GA,vy_GA,vz_GA,dens_GB,press_GB,vx_GB,vy_GB,vz_GB,&
                            dens_new,vx_new,vy_new,vz_new,press_new)
!!$                       ! BEGIN DEBUG
!!$                       if ( press_GB(i,j,k,lb) > 10000.1 ) then
!!$                          print *, 'press too low rgfm 181' , press_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( dens_GB(i,j,k,lb) - 0.125  < .000 ) then
!!$                          print *, 'dens too low rgfm 184' , dens_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vx_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vx too low rgfm 187' , vx_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vy_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vy too low rgfm 190' , vy_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vz_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vz too low rgfm 193' , vz_GB(i,j,k,lb)
!!$                       endif
!!$                       ! END DEBUG
                       ! print *, dens_new(i,j,k,lb),vx_new(i,j,k,lb),vy_new(i,j,k,lb)&
                       !      ,vz_new(i,j,k,lb),press_new(i,j,k,lb),'exiting setUpData in RGFM flB'
                    endif !if(unk(11,i,j,k,lb)<=0) then


                 endif !if ( mark(i,j,k,lb) =1)
              endif
              

              end do !k
           end do !j
        end do !i
     endif !nodetype.eq.1
  end do !lb












!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then
        do k=1+nguard,nzb+nguard
           do j=1+nguard,nyb+nguard
              do i=1+nguard,nxb+nguard
                 if (mark(i,j,k,lb) == 1) then !working in nodes immediately next to 0 level set
                    !////////////
                    !overwrite properties with properties obtained from Riemann problem (exact_MMRS_st outputs)
                    !recall "GA" vals only defined in the real region -6dx from phi =0 so far,
                    !not extrapolated into ghost A region yet
                    if (press_new(i,j,k,lb)<0.0) then
                       print *, 'press_new(i,j,k,lb) < 0 ln250 rgfm'
                    endif
                    if (unk(11,i,j,k,lb) .le. 0.0) then
                       !fluid A region; overwrite real fluid A node immediately bordering interface with star value
                       !density_new is density star


                       !debug
!!$                       if ( abs(  abs(dens_GA(i,j,k,lb))-abs(dens_new(i,j,k,lb))  )  > 0.00000000001 ) then
!!$                          print *, 'abs(dens_GA(i,j,k,lb))-abs(dens_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          print *, 'dens_GA(i,j,k,lb)',dens_GA(i,j,k,lb),'dens_new(i,j,k,lb)',dens_new(i,j,k,lb)
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(vx_GA(i,j,k,lb))-abs(vx_new(i,j,k,lb))  ) > 0.00000000001 ) then
!!$                          print *, 'abs(vx_GA(i,j,k,lb))-abs(vx_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(vy_GA(i,j,k,lb))-abs(vy_new(i,j,k,lb))  ) > 0.000000000010 ) then
!!$                          print *, 'abs(vx_GA(i,j,k,lb))-abs(vx_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(vz_GA(i,j,k,lb))-abs(vz_new(i,j,k,lb))  ) > 0.00000000001 ) then
!!$                          print *, 'abs(vz_GA(i,j,k,lb))-abs(vz_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(press_GA(i,j,k,lb))-abs(press_new(i,j,k,lb))  ) > 0.00000000001 ) then
!!$                          print *, 'abs(press_GA(i,j,k,lb))-abs(press_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
                       ! end debug
!!$                       if ( dens_GA(i,j,k,lb) .eq. dens_new(i,j,k,lb) ) then
!!$                          print *, 'ln 349 RGFM4_3Dv2 dens_GA(i,j,k,lb) .eq. dens_new(i,j,k,lb)'
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$                       endif


                       dens_GA(i,j,k,lb)  = dens_new(i,j,k,lb)
                       vx_GA(i,j,k,lb)    = vx_new(i,j,k,lb)
                       vy_GA(i,j,k,lb)    = vy_new(i,j,k,lb)
                       vz_GA(i,j,k,lb)    = vz_new(i,j,k,lb)
                       press_GA(i,j,k,lb) = press_new(i,j,k,lb)

!!$                       ! begin debug
!!$                       if ( abs(vx_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vxGA too low rgfm 278 post riemann' , vx_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vy_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vyGA too low rgfm 281 post riemann' , vy_GA(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vz_GA(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vzGA too low rgfm 284 post riemann' , vz_GA(i,j,k,lb)
!!$                       endif
!!$                       ! end debug

!!$                       ! BEGIN DEBUG
!!$                       if ( press_GA(i,j,k,lb) > 20000. ) then
!!$                          print *, 'press too high rgfm 201' , press_GA(i,j,k,lb)
!!$                       endif
!!$                       ! END DEBUG

                    end if ! (unk(11,i,j,k,lb) < 0.0)
                    !recall "GB" vals only defined in the real region +6dx from phi =0 
                    if (unk(11,i,j,k,lb).ge. 0.0) then !fluid B region
                       !fluid B region; overwrite real fluid B node immediately bordering interface with star value
                       !density_new is density star

                       !debug
!!$                       if (abs(  abs(dens_GB(i,j,k,lb))-abs(dens_new(i,j,k,lb))  )> 0.00000000001 ) then
!!$                          print *, 'abs(dens_GB(i,j,k,lb))-abs(dens_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          print *, 'dens_GB(i,j,k,lb)',dens_GB(i,j,k,lb),'dens_new(i,j,k,lb)',dens_new(i,j,k,lb)
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if (abs(  abs(vx_GB(i,j,k,lb))-abs(vx_new(i,j,k,lb))  ) >0.00000000001 ) then
!!$                          print *, 'abs(vx_GB(i,j,k,lb))-abs(vx_new(i,j,k,lb))  )> 0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(vy_GB(i,j,k,lb))-abs(vy_new(i,j,k,lb))  ) >0.00000000001 ) then
!!$                          print *, 'abs(vx_GB(i,j,k,lb))-abs(vx_new(i,j,k,lb))  )>0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(vz_GB(i,j,k,lb))-abs(vz_new(i,j,k,lb))  ) >0.00000000001 ) then
!!$                          print *, 'abs(vz_GB(i,j,k,lb))-abs(vz_new(i,j,k,lb))  ) >0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
!!$                       if ( abs(  abs(press_GB(i,j,k,lb))-abs(press_new(i,j,k,lb))  ) > 0.00000000001 ) then
!!$                          print *, 'abs(press_GB(i,j,k,lb))-abs(press_new(i,j,k,lb))  ) > 0.00000000001 '
!!$                          Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$                       endif
                       ! end debug

                       dens_GB(i,j,k,lb)  = dens_new(i,j,k,lb)
                       vx_GB(i,j,k,lb)    = vx_new(i,j,k,lb)
                       vy_GB(i,j,k,lb)    = vy_new(i,j,k,lb)
                       vz_GB(i,j,k,lb)    = vz_new(i,j,k,lb)
                       press_GB(i,j,k,lb) = press_new(i,j,k,lb)

!!$                       ! begin debug
!!$                       if ( abs(vx_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'mag vxGB too high rgfm 309 post riemann' , vx_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vy_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, ' mag vyGB too high rgfm 313 post riemann' , vy_GB(i,j,k,lb)
!!$                       endif
!!$                       if ( abs(vz_GB(i,j,k,lb)) -.01 .ge. 0. ) then
!!$                          print *, 'vzGB mag too high rgfm 315 post riemann' , vz_GB(i,j,k,lb)
!!$                       endif
!!$                       !end debug

                    endif ! (unk(11,i,j,k,lb).ge. 0.0)
                    !////////////
                 endif !(mark(i,j,k,lb)==1)
              end do !i loop
           end do !j loop
        end do !k loop

        !now overwrite the temp u values with GA, GB values
        !which are unchanged except for the values +/- 1 node from interface, and are 0 everywhere beyond their
        ! own real region and 6dx across the interface in i,j,k

     end if !nodetype.eq.1
  end do !lb


!!$
!!$  !     call extrapolation_rgfm(mype,n1,ny)
!!$  !     call extrapolation_rgfm_rev(mype,n1,ny,mark)
!!$

  !/////debug nx


  ! REPLACE unk values in band +/- 6dx so we can use PARAMESH guardcell filling
  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then
        do i=1+nguard,nxb+nguard
           do j=1+nguard,nyb+nguard
              do k=1+nguard,nzb+nguard
                 if (abs(unk(11,i,j,k,lb)).le. 4.*dx_min) then
                    ! if  (unk(11,i,j,k,lb)<0.) then
                    !if (mark(i,j,k,lb).eq.1) then
                    unk(1,i,j,k,lb) = dens_GA(i,j,k,lb)
                    unk(2,i,j,k,lb) = vx_GA(i,j,k,lb)
                    unk(3,i,j,k,lb) = vy_GA(i,j,k,lb)
                    unk(4,i,j,k,lb) = vz_GA(i,j,k,lb)
                    unk(5,i,j,k,lb) = press_GA(i,j,k,lb)
                    !endif
                    unk(6,i,j,k,lb) = dens_GB(i,j,k,lb)
                    unk(7,i,j,k,lb) = vx_GB(i,j,k,lb)
                    unk(8,i,j,k,lb) = vy_GB(i,j,k,lb)
                    unk(9,i,j,k,lb) = vz_GB(i,j,k,lb)
                    unk(10,i,j,k,lb) = press_GB(i,j,k,lb)
                    ! endif
                    ! if  (unk(11,i,j,k,lb)>0.) then
                    !if (mark(i,j,k,lb).eq.1)  then



                    !endif

                    !endif !if (unk(11,i,j,k,lb)>0.)
                 endif


              enddo
           enddo
        enddo
     endif
  enddo



  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  iopt=1  
  nlayers = nguard

  call amr_guardcell(mype,iopt,nlayers)

!!$  ! ///**  DEBUG BEGIN **/// {
!!$  do lb = 1,lnblocks
!!$     if(nodetype(lb).eq.1 .or. advance_all_levels) then
!!$        do k=1,nzb+6
!!$           do j=1,nyb+6
!!$              do i=1,nxb+6
!!$                 dx = bsize(1,lb)/real(nxb)
!!$                    dy = bsize(2,lb)/real(nyb)
!!$                    dz = bsize(3,lb)/real(nzb)
!!$                    xi = bnd_box(1,1,lb) + dx*(real(i-nguard)-.5)
!!$                    yi = bnd_box(1,2,lb) + dy*(real(j-nguard)-.5)        
!!$                    zi = bnd_box(1,3,lb) + dz*(real(k-nguard)-0.5)
!!$                 if ( abs(unk(11,i,j,k,lb)) > 2.0*dx_min ) then
!!$                 if ( abs(unk(2,i,j,k,lb))-.001 .ge. 0.) then
!!$                    print *, 'problem with x velocity RGFM before extrapolation ln 384'
!!$                    print *, 'x,y,z = ',xi,yi,zi
!!$                    print *, 'phi = ', unk(11,i,j,k,lb)
!!$                    
!!$                 endif
!!$                 if ( abs(unk(3,i,j,k,lb))-.001 .ge. 0.) then
!!$                    print *, 'problem with y velocity RGFM before extrapolation ln 387'
!!$                    print *, 'x,y,z = ',xi,yi,zi
!!$                    print *, 'phi = ', unk(11,i,j,k,lb)
!!$                 endif
!!$                 if ( abs(unk(4,i,j,k,lb))-.001 .ge. 0.) then
!!$                    print *, 'problem with z velocity RGFM before extrapolation ln 390'
!!$                    print *, 'x,y,z = ',xi,yi,zi
!!$                    print *, 'phi = ', unk(11,i,j,k,lb)
!!$                 endif
!!$                 endif
!!$              enddo
!!$           enddo
!!$        enddo
!!$     endif
!!$  enddo
!!$  ! ****** } END DEBUG


  call extrapolation_rgfmEitherSides(mype,nx,ny,nz,mark) !move properties across zero level set into ghost regions

!!$  ! ****** DEBUG
!!$  do lb = 1,lnblocks
!!$     if(nodetype(lb).eq.1 .or. advance_all_levels) then
!!$        do k=1,nzb+6
!!$           do j=1,nyb+6
!!$              do i=1,nxb+6
!!$                 dx = bsize(1,lb)/real(nxb)
!!$        dy = bsize(2,lb)/real(nyb)
!!$        dz = bsize(3,lb)/real(nzb)
!!$                 if ( abs(unk(11,i,j,k,lb)) >3.0*dx_min) then
!!$                    if ( abs(unk(2,i,j,k,lb))-.001 .ge. 0.) then
!!$                       print *, 'problem with x velocity RGFM after extrapolation'
!!$                        print *, 'ending program RGFM ln 368 with MPI_ABORT'
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$                    endif
!!$                    if ( abs(unk(3,i,j,k,lb))-.001 .ge. 0.) then
!!$                       print *, 'problem with y velocity RGFM after extrapolation'
!!$                        print *, 'ending program RGFM ln 368 with MPI_ABORT'
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$                    endif
!!$                    if ( abs(unk(4,i,j,k,lb))-.001 .ge. 0.) then
!!$                       print *, 'problem with z velocity RGFM after extrapolation'
!!$                        print *, 'ending program RGFM ln 368 with MPI_ABORT'
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$                    endif
!!$                 endif
!!$
!!$              enddo
!!$           enddo
!!$        enddo
!!$
!!$
!!$
!!$     endif
!!$  enddo
!!$
!!$  ! ****** END DEBUG
  iopt=1  
  nlayers = nguard

  call amr_guardcell(mype,iopt,nlayers)


  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then

        dx = bsize(1,lb)/real(nxb)
        dy = bsize(2,lb)/real(nyb)
        dz = bsize(3,lb)/real(nzb)
        xi = bnd_box(1,1,lb) + dx*(real(i-3)-.5)
        yi = bnd_box(1,2,lb) + dy*(real(j-3)-.5)
        zi = bnd_box(1,3,lb) + dz*(real(k-3)-.5)

        !assign now extrapolated star values in unk to dens_GA etc.
        do k=1,nzb+2*nguard
           do j=1,nyb+2*nguard
              do i=1,nxb+2*nguard
!!$                 if (abs(unk(11,i,j,k,lb))>2.0*dx_min) then
!!$                    if ( abs(unk(2,i,j,k,lb))-.001 .ge. 0.) then
!!$                       print *, 'problem with x velocity RGFM after extrapolation'
!!$                       print *, 'x,y,z = ',xi,yi,zi
!!$                    endif
!!$                    if ( abs(unk(3,i,j,k,lb))-.001 .ge. 0.) then
!!$                       print *, 'problem with y velocity RGFM after extrapolation'
!!$                       print *, 'x,y,z = ',xi,yi,zi
!!$                    endif
!!$                    if ( abs(unk(4,i,j,k,lb))-.001 .ge. 0.) then
!!$                       print *, 'problem with z velocity RGFM after extrapolation'
!!$                       print *, 'x,y,z = ',xi,yi,zi
!!$                    endif
!!$                 endif

                 if (abs(unk(11,i,j,k,lb))<= 4.0*dx_min) then
                    if (unk(11,i,j,k,lb) > 0.0 ) then ! populate ghost region of A
                       dens_GA(i,j,k,lb) = unk(1,i,j,k,lb)
                       vx_GA(i,j,k,lb) = unk(2,i,j,k,lb)
                       vy_GA(i,j,k,lb) = unk(3,i,j,k,lb)
                       vz_GA(i,j,k,lb) = unk(4,i,j,k,lb)
                       press_GA(i,j,k,lb) = unk(5,i,j,k,lb)
                    else !populate ghost region of B

                       dens_GB(i,j,k,lb) = unk(6,i,j,k,lb)
                       vx_GB(i,j,k,lb) = unk(7,i,j,k,lb)
                       vy_GB(i,j,k,lb) = unk(8,i,j,k,lb)
                       vz_GB(i,j,k,lb) = unk(9,i,j,k,lb)
                       press_GB(i,j,k,lb) = unk(10,i,j,k,lb)
                    endif
                 endif

              end do
           end do
        end do



     endif !nodetype
  enddo !lb


  unk(:,:,:,:,1:lnblocks)=t_unk(:,:,:,:,1:lnblocks) !restore original unk values
  !unk = t_unk




  do lb = 1,lnblocks
     if(nodetype(lb).eq.1 .or. advance_all_levels) then
        do k = 1+nguard,nzb+nguard
           do i=1+nguard,nxb+nguard
              do j=1+nguard,nyb+nguard
                 if ( abs(unk(11,i,j,k,lb)) .lt. (4.*dx_min) ) then
                    if (1 == -1) then  ! then run vof on interfacial cells
                       if (mark(i,j,k,lb) == 1) then
                          w2 = 1.
                          w2 = 0.
                          w1 = vof(i,j,k,lb) !working in node immediately next to interface
                          !w1 is < 0.5 by convention.
                          if (vof(i,j,k,lb) > 0.00000001) then
                             w2 = 1.0-w1 
                             if (unk(11,i,j,k,lb)<=0.0) then

                                dens_GA(i,j,k,lb) =dens_oA(i,j,k,lb)*w2+dens_GA(i,j,k,lb)*w1
                                press_GA(i,j,k,lb)=press_oA(i,j,k,lb)*w2+press_GA(i,j,k,lb)*w1
                                !note vn_GA is original values of v dot n in region A

                                !replace the smaller fraction of vdotn here with original value
                                vn_A = vx_GA(i,j,k,lb)*nx(i,j,k,lb)+vy_GA(i,j,k,lb)*ny(i,j,k,lb)+vz_GA(i,j,k,lb)*nz(i,j,k,lb)
                                vn_GA(i,j,k,lb) = (t_unk(2,i,j,k,lb)*nx(i,j,k,lb)+  t_unk(3,i,j,k,lb)*ny(i,j,k,lb)+ &
                                     t_unk(4,i,j,k,lb)*nz(i,j,k,lb))/t_unk(1,i,j,k,lb)
                                !vx_GA here has the star normal from the Riemann problem.  Replace this with the larger...
                                ! ...volume fraction (w2) of old normal (vn_GA)
                                vx_GA(i,j,k,lb)=vx_GA(i,j,k,lb)+nx(i,j,k,lb)*w2*(vn_GA(i,j,k,lb)-vn_A)      
                                vy_GA(i,j,k,lb)=vy_GA(i,j,k,lb)+ny(i,j,k,lb)*w2*(vn_GA(i,j,k,lb)-vn_A)
                                vz_GA(i,j,k,lb)=vz_GA(i,j,k,lb)+nz(i,j,k,lb)*w2*(vn_GA(i,j,k,lb)-vn_A)
                                unk(1,i,j,k,lb) = dens_GA(i,j,k,lb)
                                unk(2,i,j,k,lb) = dens_GA(i,j,k,lb)*vx_GA(i,j,k,lb)
                                unk(3,i,j,k,lb) = dens_GA(i,j,k,lb)*vy_GA(i,j,k,lb)
                                unk(4,i,j,k,lb) = dens_GA(i,j,k,lb)*vz_GA(i,j,k,lb)

                                E_cell = 0.5*(vx_GA(i,j,k,lb)**2+vy_GA(i,j,k,lb)**2+vz_GA(i,j,k,lb)**2 ) + (press_GA(i,j,k,lb) + &
                                     gamm_A*p_inf_A)/(dens_GA(i,j,k,lb)*(gamm_A-1.))
                                unk(5,i,j,k,lb) = dens_GA(i,j,k,lb)*E_cell


                             endif





                             if (unk(11,i,j,k,lb)>=0.0) then  !unk(11.. > 0
                                dens_GB(i,j,k,lb) =dens_oB(i,j,k,lb)*w2+dens_GB(i,j,k,lb)*w1
                                press_GB(i,j,k,lb)=press_oB(i,j,k,lb)*w2+press_GB(i,j,k,lb)*w1

                                vn_B = vx_GB(i,j,k,lb)*nx(i,j,k,lb)+vy_GB(i,j,k,lb)*ny(i,j,k,lb)+vz_GB(i,j,k,lb)*nz(i,j,k,lb)
                                vn_GB(i,j,k,lb) = (t_unk(7,i,j,k,lb)*nx(i,j,k,lb)+  t_unk(8,i,j,k,lb)*ny(i,j,k,lb)+ &
                                     t_unk(9,i,j,k,lb)*nz(i,j,k,lb))/t_unk(6,i,j,k,lb)

                                vx_GB(i,j,k,lb)=vx_GB(i,j,k,lb)+nx(i,j,k,lb)*w2*(vn_GB(i,j,k,lb)-vn_B)      
                                vy_GB(i,j,k,lb)=vy_GB(i,j,k,lb)+ny(i,j,k,lb)*w2*(vn_GB(i,j,k,lb)-vn_B)
                                vz_GB(i,j,k,lb)=vz_GB(i,j,k,lb)+nz(i,j,k,lb)*w2*(vn_GB(i,j,k,lb)-vn_B) 
                                unk(6,i,j,k,lb) = dens_GB(i,j,k,lb)
                                unk(7,i,j,k,lb) = dens_GB(i,j,k,lb)*vx_GB(i,j,k,lb)
                                unk(8,i,j,k,lb) = dens_GB(i,j,k,lb)*vy_GB(i,j,k,lb)
                                unk(9,i,j,k,lb) = dens_GB(i,j,k,lb)*vz_GB(i,j,k,lb)
                                E_cell = 0.5*(vx_GB(i,j,k,lb)**2+vy_GB(i,j,k,lb)**2 + vz_GB(i,j,k,lb)**2 )+(press_GB(i,j,k,lb) + &
                                     gamm_B*p_inf_B)/(dens_GB(i,j,k,lb)*(gamm_B-1.))
                                unk(10,i,j,k,lb) = dens_GB(i,j,k,lb)*E_cell

                             endif !unk(11,i,j,k,lb)>=0.
                          endif !if vof > 0.0000001


                       endif ! if mark == 1
                    endif !vof option flag if 1==1



                    if ( unk(11,i,j,k,lb) >= 0. ) then !fill in ghost

                       unk(1,i,j,k,lb) = dens_GA(i,j,k,lb)
                       unk(2,i,j,k,lb) = dens_GA(i,j,k,lb)*vx_GA(i,j,k,lb)
                       unk(3,i,j,k,lb) = dens_GA(i,j,k,lb)*vy_GA(i,j,k,lb)
                       unk(4,i,j,k,lb) = dens_GA(i,j,k,lb)*vz_GA(i,j,k,lb)

                       E_cell = 0.5*(vx_GA(i,j,k,lb)**2+vy_GA(i,j,k,lb)**2+vz_GA(i,j,k,lb)**2 ) + (press_GA(i,j,k,lb) + &
                            gamm_A*p_inf_A)/(dens_GA(i,j,k,lb)*(gamm_A-1.))
                       unk(5,i,j,k,lb) = dens_GA(i,j,k,lb)*E_cell
                    endif

                    if (unk(11,i,j,k,lb) <= 0. ) then




                       unk(6,i,j,k,lb) = dens_GB(i,j,k,lb)
                       unk(7,i,j,k,lb) = dens_GB(i,j,k,lb)*vx_GB(i,j,k,lb)
                       unk(8,i,j,k,lb) = dens_GB(i,j,k,lb)*vy_GB(i,j,k,lb)
                       unk(9,i,j,k,lb) = dens_GB(i,j,k,lb)*vz_GB(i,j,k,lb)
                       E_cell = 0.5*(vx_GB(i,j,k,lb)**2+vy_GB(i,j,k,lb)**2 + vz_GB(i,j,k,lb)**2 )+(press_GB(i,j,k,lb) + &
                            gamm_B*p_inf_B)/(dens_GB(i,j,k,lb)*(gamm_B-1.))
                       unk(10,i,j,k,lb) = dens_GB(i,j,k,lb)*E_cell
!!$                    xi = bnd_box(1,1,lb) + dx*(real(i-3)-.5)
!!$                    yi = bnd_box(1,2,lb) + dy*(real(j-3)-.5)
!!$                    zi = bnd_box(1,3,lb) + dz*(real(k-3)-.5)
!!$                    if (zi>1.9) then
!!$
!!$                       if (abs(vz_GB(i,j,k,lb)).gt.(.1)) then
!!$                          print *, 'oh no RGFM 458'
!!$                          print *, 'xyz = ',xi,yi,zi,'vz_GB = ',vz_GB(i,j,k,lb)
!!$                       endif
!!$                    endif
!!$!!!!////////////////// DEBUG
!!$                    dx = bsize(1,lb)/real(nxb)
!!$                    dy = bsize(2,lb)/real(nyb)
!!$                    dz = bsize(3,lb)/real(nzb)
!!$                    xi = bnd_box(1,1,lb) + dx*(real(i-3)-.5)
!!$                    yi = bnd_box(1,2,lb) + dy*(real(j-3)-.5)
!!$                    zi = bnd_box(1,3,lb) + dz*(real(k-3)-.5)
!!$                    if ((zi>1.9).or.(zi<.1)) then
!!$
!!$                       if ( abs( ( abs(unk(6,i,j,k,lb))-0.125 ) ).gt.0.0001 ) then
!!$                          print *, 'oh no RGFM  density before GB update 496'
!!$                          print *, 'xyz = ',xi,yi,zi,'dens = ',unk(6,i,j,k,lb)
!!$                       endif
!!$                       if ( abs( ( abs(press_GB(i,j,k,lb))-10000. ) ).gt.1. ) then
!!$                          print *, 'oh no RGFM press before GB update 500'
!!$                          print *, 'xyz = ',xi,yi,zi,'press = ',press_GB(i,j,k,lb)
!!$                       endif
!!$
!!$                    endif
!!$
!!$!!!!/////////////// END DEBUG
                    endif !unk(11 <= 0.

                 endif !abs(unk(11,.. > ...

!!$                 if( abs(unk(11,i,j,k,lb))< 6.0*dx_min) then
!!$                    if (unk(11,i,j,k,lb) < 0.0) then
!!$                       !        vnorm(i,j,k,lb) = vx_GA(i,j,k,lb)*nx(i,j,k,lb)+vy_GA(i,j,k,lb)*ny(i,j,k,lb)
!!$                    else
!!$                       !        vnorm(i,j,k,lb) = vx_GB(i,j,k,lb)*nx(i,j,k,lb)+vy_GB(i,j,k,lb)*ny(i,j,k,lb)
!!$                    endif
!!$                 endif
!!$
!!$
!!$                 if (abs(unk(11,i,j,k,lb)) < 2.0*dx_min) then
!!$                    !       print *,vx_GA(i,j,k,lb),vx_GB(i,j,k,lb),vnorm(i,j,k,lb),nx(i,j,k,lb)
!!$                 endif



              end do !j
           end do !i
        enddo !k

     endif !nodetype
  enddo !lb

  !    call extrapolation_vnorm(mype,nx,ny,vnorm,mark)   
  !    if(mod(steps,20).eq.0) call check_mass(mype,mark,vof,steps)





  return
end subroutine RGFM




subroutine setUpDataForRiemannProblem(xii,yii,zii,fluid,i,j,k,lb,nx,ny,nz,curvature,&
     dens_GA,press_GA,vx_GA,vy_GA,vz_GA,dens_GB,press_GB,vx_GB,vy_GB,vz_GB,&
     dens_new,vx_new,vy_new,vz_new,press_new)

  use paramesh_dimensions
  use physicaldata
  use tree
  include 'mpif.h'


  real :: dens_IL,vn_I,press_IL,dens_IR,press_IR,&
       xI,yI,zI,&
       dens_A,vx_A,vy_A,vz_A,vn_A,vt_A,press_A, &
       dens_B,vx_B,vy_B,vz_B,vn_B,vt_B,press_B,E_cell, &
       gavgA_A,gavgB_B,Lavg_A,Lavg_B,tem_A,tem_B
  real :: p,q,r,dx,dy,dz,w1,w2,qCond,nxA
  real :: kappa,mvapFlux,molFrac_B2
  integer :: ii,jj

  character(len=1), intent(in) :: fluid

  integer,intent(in) :: i,j,k,lb
  real,intent(in) ::  nx(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       ny(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       nz(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       curvature(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)

  real,intent(in) :: dens_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       press_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vx_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),& 
       vy_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vz_GA(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks), &
       dens_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       press_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vx_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vy_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks),&
       vz_GB(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)

  real,intent(in) :: xii,yii,zii

  real, intent(out) :: dens_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)
  real, intent(out) :: press_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)
  real, intent(out) :: vx_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)
  real, intent(out) :: vy_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)
  real, intent(out) :: vz_new(nxb+2*nguard,nyb+2*nguard,nzb+2*nguard,maxblocks)
  real::phiVal
  real :: normalizedDistance





  !this function uses the values in a band +/- 2*nguarddx away from the 0 level set in all directions,
  !to ultimately return the star region values of the Reimann problem at the interface,
  !and assign them to the neighbor nodes in A or B.
  !this function should only be called at the 1st nodes immediately away from the interface.
  phiVal=unk(11,i,j,k,lb)
  !step1: move along normal direction to find zero level set (gives xI fractional node #,  eg node xI = 1.2 )
  xI = real(i)- unk(11,i,j,k,lb)*nx(i,j,k,lb)/dx_min



  yI = real(j)  - unk(11,i,j,k,lb)*ny(i,j,k,lb)/dx_min
  zI = real(k)  - unk(11,i,j,k,lb)*nz(i,j,k,lb)/dx_min
  normalizedDistance = unk(11,i,j,k,lb)/dx_min
  !step2A: now move from this 0 level set coordinate 1.5 units in -1*nhat dir. These coordinates are p,q,r
  p = xI - 1.733*nx(i,j,k,lb) !1.75 = sqrt(3)+epsilon. This guarantees that the convex hull for interpolation...
  !...does not cross interface.




  q = yI - 1.733*ny(i,j,k,lb)
  r = zI - 1.733*nz(i,j,k,lb)
  if (isnan(p).or.isnan(q).or.isnan(r)) then
     print *,'nan ln 928 RGFM4'
     
     print *, 'p = ',p,'q= ',q,'r = ',r
     print *, 'p = ',p,'q = ',q,'r = ',r
     print *, 'nxA = ', nx(i,j,k,lb), 'nyA =', ny(i,j,k,lb), 'nzA =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  !DEBUG{
  if ((p.le.1) .or. (p.gt.10).or.(q.le.1) .or. (q.gt.10).or.(r.le.1) .or. (r.gt.10) ) then
     print *,'RGFM out of bounds'    
     print *, 'p = ',p,'q = ',q,'r = ',r
     print *, 'nxA = ', nx(i,j,k,lb), 'nyA =', ny(i,j,k,lb), 'nzA =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  !}END DEBUG
  !populate L of discontinuity properties at the 1.5dx points using trilinear interpolation         
  dens_A = trilinear_interpolation(p,q,r,dens_GA,lb,2)
  press_A= trilinear_interpolation(p,q,r,press_GA,lb,2)
  vx_A   = trilinear_interpolation(p,q,r,vx_GA,lb,0)
  vy_A   = trilinear_interpolation(p,q,r,vy_GA,lb,0)
  vz_A   = trilinear_interpolation(p,q,r,vz_GA,lb,0)

  !debug{
!!$  if ( (  abs(dens_GA(i,j,k,lb))-abs(dens_A)  ) .ne. 0.000000000000 ) then
!!$     print *, 'abs(dens_GA(i,j,k,lb))-abs(dens_A(i,j,k,lb))  ) .ne. 0.000000000000 '
!!$     print *, 'dens_GA(i,j,k,lb)',dens_GA(i,j,k,lb),'dens_A',dens_A
!!$     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$  endif
!!$  if ( (  abs(dens_A-1.0000000000000)  ) >0.000000000001 ) then
!!$     print *, 'abs(dens_A-1.0000000000000)  ) >0.000000000001 ) '
!!$     print *, 'dens_GA(i,j,k,lb)',dens_GA(i,j,k,lb),'dens_A',dens_A
!!$     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$
!!$  endif

  if (dens_A == 3.14159265359) then
     print *, 'RGFM called from fluid ',fluid
     print *, 'error occurred in A values'
     print *, 'nxA = ', nx(i,j,k,lb), 'nyA =', ny(i,j,k,lb), 'nzA =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  if (isnan(dens_A)) then
     print *, 'dens_A isnan rgfm ln 991'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(press_A)) then
     print *, 'press_A isnan rgfm ln 991'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  
  if (press_A == 23.14159265359) then
     print *, 'press_A problem'
     print *, 'RGFM called from fluid ',fluid
     print *, 'error occurred in A values'
     print *, 'nxA = ', nx(i,j,k,lb), 'nyA =', ny(i,j,k,lb), 'nzA =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  !} end debug


  !calculate v dot n for region A 1.5dx away from phi = 0 (v_L in the reimann problem)
  vn_A   = vx_A*nx(i,j,k,lb)+vy_A*ny(i,j,k,lb)+vz_A*nz(i,j,k,lb)
!!$  if (abs(vn_A) > 0.1 ) then
!!$  print *, 'vn_A mag too big',vn_A
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$  endif

  !step2B: now move from this 0 level set coordinate 1.5 units in -1*nhat dir. These coordinates are p,q,r
  p = xI + 1.733*nx(i,j,k,lb)
!!$  !if (p.le.0.0) print *, 'nxB = ', nx(i,j,k,lb)
!!$  !if (p.gt.100.) print *, 'nxB = ', nx(i,j,k,lb), 'nyB =', ny(i,j,k,lb), 'nzB =', nz(i,j,k,lb),&
!!$   !    'x= ', xii, 'y= ', yii, 'z= ', zii
  q = yI + 1.733*ny(i,j,k,lb)
  r = zI + 1.733*nz(i,j,k,lb)
  !DEBUG{
   if (isnan(p).or.isnan(q).or.isnan(r)) then
     print *,'nan ln 928 RGFM4'
     
     print *, 'p = ',p,'q= ',q,'r = ',r
     print *, 'p = ',p,'q = ',q,'r = ',r
     print *, 'nxA = ', nx(i,j,k,lb), 'nyA =', ny(i,j,k,lb), 'nzA =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  if ((p.le.1) .or. (p.gt.10).or.(q.le.1) .or. (q.gt.10).or.(r.le.1) .or. (r.gt.10) ) then
     print *, 'p = ',p,'q = ',q,'r = ',r

     print *, 'nxA = ', nx(i,j,k,lb), 'nyA =', ny(i,j,k,lb), 'nzA =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  !}END DEBUG


  !populate R of discontinuity properties at the 1.5dx points using trilinear interpolation            
  dens_B = trilinear_interpolation(p,q,r,dens_GB,lb,2)

  press_B= trilinear_interpolation(p,q,r,press_GB,lb,2)
  
  vx_B   = trilinear_interpolation(p,q,r,vx_GB,lb,0)
  vy_B   = trilinear_interpolation(p,q,r,vy_GB,lb,0)
  vz_B = trilinear_interpolation(p,q,r,vz_GB,lb,0)
  !debug{
!!$   if ( (  abs(dens_B-1.0000000000000)  ) >0.000000000001 ) then
!!$     print *, 'abs(dens_B-1.0000000000000)  ) >0.000000000001 ) '
!!$     print *, 'dens_GA(i,j,k,lb)',dens_GB(i,j,k,lb),'dens_B',dens_B
!!$     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$  endif
  !}end debug
  if (dens_B == 3.14159265359) then
     print *, 'dens_B problem RGFM called from fluid ',fluid
     print *, 'error occurred in B values'
     print *, 'nxB = ', nx(i,j,k,lb), 'nyB =', ny(i,j,k,lb), 'nzB =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif
  if (press_B == 23.14159265359) then
     print *, 'press_B problem  RGFM called from fluid ',fluid
     print *, 'RGFM called from fluid ',fluid
     print *, 'error occurred in B values'
     print *, 'nxB = ', nx(i,j,k,lb), 'nyB =', ny(i,j,k,lb), 'nzB =', nz(i,j,k,lb),&
          'xI= ', xI, 'yI= ', yI, 'z= ', zI , 'i,j,k= ',i,j,k,'phi = ', phiVal
     print *, 'phi/dx = ',normalizedDistance,'x= ', xii, 'y= ', yii, 'z= ', zii
  endif


  !curvature_B   = trilinear_interpolation(p,q,r,curvature,lb)

  !calculate v dot n for region B 1.5dx away from phi = 0  (v_R in the riemann problem)
  vn_B   = vx_B*nx(i,j,k,lb)+vy_B*ny(i,j,k,lb)+vz_B*nz(i,j,k,lb)

!!$  if (abs(vn_B) > 0.1 ) then
!!$  print *, 'vn_B mag too big',vn_B
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$  endif



  !         call exact_MMRS(dens_A,vn_A,press_A,dens_B,vn_B,press_B,dens_IL,vn_I, &
  !                         press_IL,dens_IR,press_IR,xii,yii)


  kappa=curvature(i,j,k,lb) !use the curvature at the zero level set itself (no need to move 1.5dx)
  ! kappa = 0.
  !"I" refers to star region values.   
  call exact_MMRS_st(j,dens_A,vn_A,press_A,gamm_A,dens_B,vn_B,press_B,&
       gamm_B,kappa,vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR,xii,yii,zii)
  
  
  

!!$  if (abs(vn_I) > 1. ) then
!!$  print *, 'vn_I mag too big',vn_I
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$  endif

  if(1==1) then             
     if(fluid=='B') then
        !assign the star values obtained from exact_MMRS.. to real interfacial fluid B node
        !immediately next to interface.

!!$        !debug code{
!!$        dens_new(i,j,k,lb)  = dens_B
!!$
!!$        ! replace normal component of velocity at nodes immediately next to  with star normal
!!$        vx_new(i,j,k,lb)    = vx_B
!!$        !+ (vn_I-vn_B)*nx(i,j,k,lb)
!!$
!!$        vy_new(i,j,k,lb)    = vy_B
!!$        !+ (vn_I-vn_B)*ny(i,j,k,lb)
!!$        vz_new(i,j,k,lb)    = vz_B
!!$        !+ (vn_I-vn_B)*nz(i,j,k,lb)
!!$        press_new(i,j,k,lb) = press_B
!!$        !end debug }


        !real code {
        !compute the normal velocity of the calling cell and overwrite the sqrt3 interp val stored in vn
        vn_B = vx_GB(i,j,k,lb)*nx(i,j,k,lb)+vy_GB(i,j,k,lb)*ny(i,j,k,lb)+vz_GB(i,j,k,lb)*nz(i,j,k,lb)
        dens_new(i,j,k,lb)  = dens_IR
     
         
        ! replace normal component of velocity at nodes immediately next to  with star normal
        vx_new(i,j,k,lb)    = vx_GB(i,j,k,lb) + (vn_I-vn_B)*nx(i,j,k,lb) 
        vy_new(i,j,k,lb)    = vy_GB(i,j,k,lb) + (vn_I-vn_B)*ny(i,j,k,lb)
        vz_new(i,j,k,lb)    = vz_GB(i,j,k,lb) + (vn_I-vn_B)*nz(i,j,k,lb)
        press_new(i,j,k,lb) = press_IR
        if (press_new(i,j,k,lb) <(0.01)) then
           print *, 'press_IR too small, fluid: ',fluid
            Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
           endif
       

!!$        !debug{
!!$        if (abs(vx_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 931'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vy_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 935'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vz_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 939'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        !enddebug}
           
        

        
!!$        if (dens_IR >100.) print *, 'dens too large, fluid: ',fluid
!!$        if (press_new(i,j,k,lb) <(0.0)) print *, 'press_IR negative, fluid: ',fluid
!!$         if (isnan(dens_IR)) print *, 'dens_IR isnan, fluid: ',fluid
!!$        if (isnan(press_IR) ) print *, 'press_IR isnan, fluid: ',fluid
     end if !if(fluid=='B') then

     if(fluid=='A') then
        !assign the star values obtained from exact_MMRS.. to real interfacial fluid A
        !real code{
        !compute the normal velocity of the calling cell and overwrite the sqrt3 interp val stored in vn
        vn_A = vx_GA(i,j,k,lb)*nx(i,j,k,lb)+vy_GA(i,j,k,lb)*ny(i,j,k,lb)+vz_GA(i,j,k,lb)*nz(i,j,k,lb)
        dens_new(i,j,k,lb)  = dens_IL
        
        vx_new(i,j,k,lb)    = vx_GA(i,j,k,lb) + (vn_I-vn_A)*nx(i,j,k,lb) 
        vy_new(i,j,k,lb)    = vy_GA(i,j,k,lb) + (vn_I-vn_A)*ny(i,j,k,lb)
        vz_new(i,j,k,lb)    = vz_GA(i,j,k,lb) + (vn_I-vn_A)*nz(i,j,k,lb)
        press_new(i,j,k,lb) = press_IL

        
        !if (dens_IL >100.) print *, 'dens too large, fluid: ',fluid
        if (press_new(i,j,k,lb) <(0.01)) then
           print *, 'press_IL too small, fluid: ',fluid
          
            Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
           endif
           
        if (isnan(dens_IL)) print *, 'dens_IL isnan, fluid: ',fluid
        if (isnan(press_IL) ) print *, 'press_IL isnan, fluid: ',fluid
         !debug{
!!$        if (abs(vx_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 968'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vy_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 972'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vz_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 976'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
        !enddebug}

        !end real code}

        
!!$        !debug code{
!!$        
!!$        dens_new(i,j,k,lb)  = dens_A
!!$
!!$        vx_new(i,j,k,lb)    = vx_A
!!$        !+ (vn_I-vn_A)*nx(i,j,k,lb) 
!!$        vy_new(i,j,k,lb)    = vy_A
!!$        !+ (vn_I-vn_A)*ny(i,j,k,lb)
!!$        vz_new(i,j,k,lb)    = vz_A
!!$        !+ (vn_I-vn_A)*nz(i,j,k,lb)
!!$        press_new(i,j,k,lb) = press_A
!!$        !end debug code }

!!$        if (abs(vx_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 1075'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vy_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 1079'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vz_new(i,j,k,lb))>0.000000000001) then
!!$           print *, 'rgfm4 1083'
!!$           Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif

     end if !if(fluid=='A') then
  endif !if(1==1) then  





end subroutine setUpDataForRiemannProblem

   
