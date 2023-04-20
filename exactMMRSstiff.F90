!subroutine exact_MMRS(jj,dens_L,vn_L,press_L,dens_R,vn_R,press_R,vn_I,dens_IL,vn_IL,&
!    press_IL,dens_IR,vn_IR,press_IR,gavgA_A,gavgB_B)

! subroutine exact_MMRS(jj,dens_L,vn_L,press_L,gavgA_A,dens_R,vn_R,press_R,gavgB_B,&
! molFrac_B2,kappa,Tii,TiiL,TiiR,mvapFlux,vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR,prtThis)

subroutine exact_MMRS_st(jj,dens_L,vn_L,press_L,gavgA_A,dens_R,vn_R,press_R,gavgB_B,&
     kappa,vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR,xii,yii,zii)


  use paramesh_dimensions
  use physicaldata
  include 'mpif.h'

  integer :: i,jj,i_index
  !    integer,intent(in) :: prtThis
  ! real,intent(in)::dens_L,vn_L,press_L,gavgA_A,dens_R,vn_R,press_R,gavgB_B,&
  ! molFrac_B2,kappa,Tii,TiiL,TiiR
  ! real,intent(out)::mvapFlux,vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR

  real,intent(in)::dens_L,vn_L,press_L,gavgA_A,dens_R,vn_R,press_R,gavgB_B,kappa,&
       xii,yii,zii
  real,intent(out)::vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR
  real :: mvapFlux

  character(len=2) :: cas
  real :: SL,ML,SR,MR,STL,SHL,STR,SHR,pI,vI,rIL,rIR,xc,rcell,vcell,pcell  
  real :: t,gamL,gamR,BL,BR,AL,AR,r0L,r0R,diap,rb,pDiff,uDiff
  real :: WL(1:3),WR(1:3)    
  real :: sigma,Temperature,T_gas,T_liq,psat,c1,C,c2,pv,rdot,udis,Pdis
  if (abs(press_R).gt.(10**13)) then
     print *, 'Pressure diverges'
  endif


  !sigma=0.5312056307952*(10.0e-6)/(4.8e-3)*10.0*1000.0/942.0*0.01 !0.0728
  sigma = 0.0728 !water air sim value
  !sigma = 800.
  mvapFlux=0.0 !evaporation

  ! Temperature=Tii; !108.0;
  ! T_gas=TiiR; !Temperature;
  ! T_liq=TiiL; !Temperature;
  ! psat=75000.0
  ! c1=(dens_R/dens_L)**(1.0/3.0);
  ! C=(1.0-c1)*exp(-0.5/(1.0/c1-1.0));
  ! c2=(2.0*C/(2.0-C))*sqrt(wB2/(2.0*3.14*8.314*1000.0));

  ! pv=molFrac_B2*press_R
  ! mvapFlux=c2*(psat/sqrt(T_liq)-pv/sqrt(T_gas)) !*1000.0
  !mvapFlux=c2*(psat/sqrt(Tii)-pv/sqrt(Tii)) !*1000.0

  if (dens_R < .000000001) then
     print *, 'dens_R <0 ln 58 exactMMRSstiff,dens_R= '
     print *, dens_R
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (dens_L < .000000001) then
     print *, 'dens_L <0 ln 58 exactMMRSstiff,dens_L= '
     print *, dens_L
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (press_R < .000000001) then
     print *, 'press_R <0 ln 58 exactMMRSstiff,press_R= '
     print *, press_R
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (press_L < .000000001) then
     print *, 'press_L <0 ln 58 exactMMRSstiff,press_L= '
     print *,press_L
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(press_R)) then
     print *, 'press_R Nan ln 58 exactMMRSstiff,press_R= '
     print *, press_R
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(press_L)) then
     print *, 'press_L Nan ln 58 exactMMRSstiff,press_L= '
     print *,press_L
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif


  
  rdot=mvapFlux*(1.0/dens_L)  
  udis=mvapFlux*(1.0/dens_R-1.0/dens_L)
  Pdis=(-sigma*kappa-mvapFlux*udis)

!!$  if ( (kappa .lt. 12.5) .or. (kappa .gt. 14.) ) then
!!$print *, '( (kappa .lt. 6.) .or. (kappa .gt. 7.) exact MMRS'
!!$print *, 'kappa = ',kappa
!!$print *, 'x,y,z = ',xii,yii,zii
!!$ Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$endif


  
  if (isnan(dens_L) ) then
     print *, 'dens_L NaN ln 58 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(dens_R) ) then
     print *, 'dens_R NaN ln 58 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(udis) ) then
     print *, 'udis NaN ln 58 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(kappa) ) then
     print *, 'kappa NaN ln 58 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(Pdis) ) then
     print *, 'Pdis NaN ln 58 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (press_R < 0.05) then
      print *, 'Pdis NaN ln 103 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  ! mvapFluxG=mvapFlux
  ! uDisGlobal=udis

  !if(jj==4) print *, 'mvapFlux,molFrac_B2,pv,T_liq,T_gas',mvapFlux,molFrac_B2,pv,T_liq,T_gas


  rb=0.0;
!!$  pDiff= -0.1
  uDiff=0.0;

!!$  rb=rdot;
  pDiff=Pdis;
!!$  uDiff=udis;

  ! rb=0.10750003558678385        ;
  ! pDiff=-12557.882821483616;
  ! uDiff=116.81747594721256       ;      

  ! AL=AAo; r0L=D0; BL=BBO;        

  AL=0.0; r0L=0.0; BL=p_inf_A;        
  AR=0.0;r0R=0.0; BR=0.0;

  gamL=gavgA_A; rL=dens_L; vL=vn_L; pL=press_L;
  WL=(/rL, vL, pL/);

  gamR=gavgB_B; rR=dens_R; vR=vn_R; pR=press_R;  
  WR=(/rR, vR, pR/);
!!$    ! BEGIN DEBUG
!!$     if (abs(vR) > 0.1 ) then
!!$        print *, 'vR mag too big exactMMRS',vR
!!$        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vL) > 0.1 ) then
!!$  print *, 'vL mag too big exactMMRS',vL
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$endif
!!$! END DEBUG

  call getStarValues(cas,pI,pIL,pIR,vI,vIL,vIR,rIL,rIR,WL,WR,gamL,gamR,BL,BR,AL,AR,r0L,r0R,rb,pDiff,uDiff,&
       xii,yii,zii)

  vn_I=vI; dens_IL=rIL; vn_IL=vIL; press_IL=pIL; dens_IR=rIR; vn_IR=vIR; press_IR=pIR;

!!$    ! BEGIN DEBUG
!!$     if (abs(vI) > 0.1 ) then
!!$        print *, 'vI mag too big exactMMRS ln 107',vI
!!$        Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$        endif
!!$        if (abs(vIL) > 0.1 ) then
!!$  print *, 'vL mag too big exactMMRS',vIL
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$endif
!!$if (abs(vIR) > 0.1 ) then
!!$  print *, 'vL mag too big exactMMRS',vIL
!!$  Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$endif
!!$
!!$! END DEBUG

if (press_R < 0.05) then
      print *, 'Pdis NaN ln 167 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(pIR) ) then
     print *, 'pIR NaN'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif

  !    vn_I=0.10750003558678385 ;
  !    dens_IL= 1.000026674548774e+03; vn_IL=-0.041412774758583; press_IL=1.637123273656890e+05;
  !    dens_IR= 1.275496313287438; vn_IR=1.164607638234140e+02; press_IR=1.512471784737372e+05;


  !print *, 'mvapFlux=',mvapFlux
  if(1==-1) then
     print *, 'mvapFlux=',mvapFlux
     print *, 'psat,T_liq,pv,T_gas=',psat,T_liq,pv,T_gas
     print *, 'rdot,udis,Pdis=',rdot,udis,Pdis
     !print *, 'gavgA_A,gavgB_B=',gavgA_A,gavgB_B
     !print *, 'vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR',vn_I,dens_IL,vn_IL,press_IL,dens_IR,vn_IR,press_IR

     print *, 'vn_I=',vn_I
     print *,'rL, vL, pL=',rL, vL, pL
     print *, 'dens_IL,vn_IL,press_IL=',dens_IL,vn_IL,press_IL
     print *,'rR, vR, pR=',rR, vR, pR
     print *, 'dens_IR,vn_IR,press_IR=',dens_IR,vn_IR,press_IR

  endif



!!$  if (abs(pIR)<1.0e20) then
!!$  else
!!$     print *, 'pIR NaN'
!!$  endif

end subroutine exact_MMRS_st









subroutine getStarValues(cas,pI,pIL,pIR,vI,vIL,vIR,rIL,rIR,WL,WR,gamL,gamR,BL,BR,AL,AR,r0L,r0R,rb,pDiff,uDiff,&
     xii,yii,zii)

  real, intent (out) :: pI,pIL,pIR,vI,vIL,vIR,rIL,rIR
  character(len=2), intent(out) :: cas
  real, intent (in) :: gamL,gamR,BL,BR,AL,AR,r0L,r0R,rb,pDiff,uDiff,xii,yii
  real, intent (in) :: WL(1:3),WR(1:3),zii

  real :: omg,TOL
  integer :: i


  omg=.51;
  TOL=1.0e-13;

  rL=WL(1); vL=WL(2); pL=WL(3);
  rR=WR(1); vR=WR(2); pR=WR(3);
!!$  if (pL .eq. pR) then
!!$     print *, 'pL=pR'
!!$     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
!!$     endif

  pI=max(.0001,0.5*(pL+pR-3.0*pDiff)+0.0);
  !pI = 100.
  !    print *, 'pI=',pI
  pIL=pI;
  pIR=pI+pDiff;
  if (isnan(pIL) ) then
     print *, 'pIL NaN ln 247 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(pDiff) ) then
     print *, 'pDiff NaN ln 251 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  if (isnan(pIR) ) then
     print *, 'pIR NaN ln 255 exactMMRSstiff'
     Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
  endif
  do i = 1,3000

     if (pIL<=pL  .and.  pIR<=pR) then
        cas='RR';
        !            print *,'pL pI pIL pIR pR=',pL, pI, pIL, pIR, pR
        store=pI;                          
        call fRarStiff(fLrar,fdLrar,pIL,pL,rL,BL,gamL);
        call fRarStiff(fRrar,fdRrar,pIR,pR,rR,BR,gamR);

        f=fLrar+fRrar+vR-vL-uDiff;
        f_prime=fdLrar+fdRrar;

        pI=pI-omg*f/f_prime; pIL=pI; pIR=pI+pDiff;
        vI=0.5*(vR+vL)+0.5*(fRrar-fLrar)-rb-0.5*uDiff; vIL=vI+rb; vIR=vIL+uDiff;
        rIL=rL*((pIL+BL)/(pL+BL))**(1/gamL);
        rIR=rR*((pIR+BR)/(pR+BR))**(1/gamR);

     elseif (pIL<=pL  .and.  pIR>pR) then
        cas='RS';
        !            print *,'pL pI pIL pIR pR=',pL, pI, pIL, pIR, pR
        store=pI;            
        call fRarStiff(fLrar,fdLrar,pIL,pL,rL,BL,gamL);
        call fShkStiff(fRshk,fdRshk,pIR,pR,rR,BR,gamR);

        f = fLrar+fRshk+vR-vL-uDiff;
        f_prime=fdLrar+fdRshk;

        pI=pI-omg*f/f_prime; pIL=pI; pIR=pI+pDiff;
        vI=0.5*(vR+vL)+0.5*(fRshk-fLrar)-rb-0.5*uDiff; vIL=vI+rb; vIR=vIL+uDiff;
        rIL=rL*((pIL+BL)/(pL+BL))**(1/gamL);
        c1=(pIR+BR)/(pR+BR);
        c2=(gamR-1)/(gamR+1);
        rIR=rR*(c1+c2)/(c1*c2+1);

     elseif (pIL>pL  .and.  pIR<=pR) then
        cas='SR';
        !            print *,'pL pI pIL pIR pR=',pL, pI, pIL, pIR, pR
        store=pI;
        call fShkStiff(fLshk,fdLshk,pIL,pL,rL,BL,gamL);
        call fRarStiff(fRrar,fdRrar,pIR,pR,rR,BR,gamR);

        f=fLshk+fRrar+vR-vL-uDiff;
        f_prime=fdLshk+fdRrar;

        pI=pI-omg*f/f_prime; pIL=pI; pIR=pI+pDiff;
        vI=0.5*(vR+vL)+0.5*(fRrar-fLshk)-rb-0.5*uDiff; vIL=vI+rb; vIR=vIL+uDiff;
        c1=(pIL+BL)/(pL+BL);
        c2=(gamL-1)/(gamL+1);
        rIL=rL*(c1+c2)/(c1*c2+1);
        rIR=rR*((pIR+BR)/(pR+BR))**(1/gamR);

     elseif (pIL>pL  .and.  pIR>pR) then
        cas='SS';
        !            print *,'pL pI pIL pIR pR=',pL, pI, pIL, pIR, pR,gamR
        store=pI;
        call fShkStiff(fLshk,fdLshk,pIL,pL,rL,BL,gamL);
        call fShkStiff(fRshk,fdRshk,pIR,pR,rR,BR,gamR);

        f=fLshk+fRshk+vR-vL-uDiff;
        f_prime=fdLshk+fdRshk;

        pI=pI-omg*f/f_prime; pIL=pI; pIR=pI+pDiff;        
        vI=0.5*(vR+vL)+0.5*(fRshk-fLshk)-rb-0.5*uDiff; vIL=vI+rb; vIR=vIL+uDiff;
        c1=(pIL+BL)/(pL+BL);
        c2=(gamL-1)/(gamL+1);
        rIL=rL*(c1+c2)/(c1*c2+1);
        c1=(pIR+BR)/(pR+BR);
        c2=(gamR-1)/(gamR+1);
        rIR=rR*(c1+c2)/(c1*c2+1);

     end if

     if(i .eq. 2999) print *, '---------Too many iterations at :',xii,yii,zii
     if(i .eq. 2999) print *,'pI = ',pI,'pIL = ',pIL,'pIR = ',pIR,'pL = ',pL,'pR = ',pR
     if(i .eq. 2999) print *,'rI = ',rI,'rIL = ',rIL,'rIR = ',rIR,'rL = ',rL,'rR = ',rR
     if(i .eq. 2999) print *,'vI = ',vI,'vIL = ',vIL,'vIR = ',vIR,'vL = ',vL,'vR = ',vR
     if (i .eq. 2999) Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
     if(2.0*abs(pI-store)/(pI+store) < TOL) goto 123
     !        print *,'cas=',cas
     !        print *,'i=',i


     if(isnan(pI) .or. isnan(pIL) .or. isnan(pIR) .or. &
     isnan(rI) .or. isnan(rIL) .or. isnan(rIR) .or. &
     isnan(vI) .or. isnan(vIL) .or. isnan(vIR) ) then
        print *, 'NaN found in exactMMRSstiff ln 343, iteration# ',i,' at :',xii,yii,zii
        print *, 'Pdiff = ',pdiff
      print *,'pI = ',pI,'pIL = ',pIL,'pIR = ',pIR,'pL = ',pL,'pR = ',pR
      print *,'rI = ',rI,'rIL = ',rIL,'rIR = ',rIR,'rL = ',rL,'rR = ',rR
      print *,'vI = ',vI,'vIL = ',vIL,'vIR = ',vIR,'vL = ',vL,'vR = ',vR
      Call MPI_ABORT(MPI_COMM_WORLD,ierrorcode,ierr)
      endif





  end do

123 continue

end subroutine getStarValues








subroutine fShkStiff(fK,fdK,pI,pK,rK,BK,gamK)

    real, intent(out) :: fK,fdK
    real, intent(in) :: pI,pK,rK,BK,gamK
    real :: pIBar,pKBar,bKo,aK
 
    pIBar=pI+BK;
    pKBar=pK+BK;
    bKo=(gamK-1.0)/(gamK+1.0)*pKBar;
    aK=2.0/(gamK+1.0)/rK;
    fK=sqrt(aK/(pIBar+bKo))*(pI-pK);
    fdK=sqrt(aK/(pIBar+bKo))*(1.0-(pI-pK)/2.0/(pIBar+bKo));

end subroutine fShkStiff


subroutine fRarStiff(fK,fdK,pI,pK,rK,BK,gamK)

    real, intent(out) :: fK,fdK
    real, intent(in) :: pI,pK,rK,BK,gamK
    real :: pIBar,pKBar,cK,trm1,exp1,exp2

    pIBar=pI+BK;
    pKBar=pK+BK;
    cK=sqrt(gamK*pKBar/rK);

    trm1=2.0*cK/(gamK-1.0);
    exp1=(gamK-1.0)/2.0/gamK;
    fK=trm1*((pIBar/pKBar)**exp1-1.0);

    exp2=-(gamK+1.0)/2.0/gamK;
    fdK=1.0/rK/cK*(pIBar/pKBar)**exp2;

end subroutine fRarStiff


subroutine StiffShockSpeed(SK,MK,pI,pK,rK,vK,BK,gamK,dir)

    real, intent(out) :: SK,MK
    real, intent(in) :: pI,pK,rK,vK,BK,gamK
    integer, intent(in) :: dir
    real :: pIBar,pKBar,bKo,aK,cK,Q

    pIBar=pI+BK;
    pKBar=pK+BK;
    bKo=(gamK-1.0)/(gamK+1.0)*pKBar;
    aK=2.0/(gamK+1.0)/rK;
    cK=sqrt(gamK*pKBar/rK);
    Q=sqrt((pIBar+bKo)/aK);
    SK=vK+dir*Q/rK;
    MK=SK/cK;

end subroutine StiffShockSpeed


subroutine StiffRarefactionSpeed(STK,SHK,cK,pI,rIK,vI,pK,rK,vK,BK,gamK,dir)

    real, intent(out) :: STK,SHK,cK
    real, intent(in) :: pI,rIK,vI,pK,rK,vK,BK,gamK
    integer, intent(in) :: dir
    real :: pIBar,cIK,pKBar

    pIBar=pI+BK;
    cIK=sqrt(gamK*pIBar/rIK);
    STK=vI+dir*cIK;

    pKBar=pK+BK;
    cK=sqrt(gamK*pKBar/rK);
    SHK=vK+dir*cK;


end subroutine StiffRarefactionSpeed



