 SUBROUTINE Add_Aggregate(pa, harvest, key, primaryProd, TEPProd, sedProd, oilProd)
  
! This creates new primary particles according to the key:
! key<=12: phytoplankton (12 pp particles)
! key<=15: TEP  (3 TEP particles)
! key<=17: sediment (2 sed part) 
! key<=20: oil  (3 oil particles)

! input is in terms of moles of organic carbon.
! n is #/m3
   USE the_info
   USE  data_sim
   USE data_stations
   
   IMPLICIT NONE
   TYPE(agg_part), INTENT(INOUT) :: pa
   INTEGER, INTENT(IN) :: key
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, primaryProd, sedProd, TEPProd, oilProd
   DOUBLE PRECISION :: tvol, age
   INTEGER :: i, phyto_n, TEP_n, sedi_n, oil_n
   
   ! diatoms contain 15pmol(ie 10^(-12)) orgC, 5pmolBSi
   ! coccolith contain 7pmol orgC, 5pmolCaCO3
   ! dust contain 1pmol kaolinite (clay)
   ! primaryProd [mgC/m2/timestep]
   ! phyto_n is the number of phytoplankton super particles produced in a timestep
 
   pa%orgC = 0.0
   pa%mineral = 0.0
   pa%TEP = 0.0
   pa%oil = 0.0
   pa%o = 0.0
   pa%sedi = 0.0
   pa%d = 0.0
   pa%s = 0.0
   pa%b = 0.0
   pa%af = 0.0
   ! Depth layer
   pa%idl = 1
   !Add a time stamp to the aggragate
   pa%t = 0.0
   pa%visc = 0.0
   

   ! Define a fractal dimension
   !pa%frac = 2.0
   CALL fractal_dimension(pa)   

   phyto_n = 12
   TEP_n = 3
   sedi_n = 2
   oil_n = 3
   IF (key .le. phyto_n) THEN
      CALL OrgC_Primary_Particle(pa, harvest, primaryProd, phyto_n)
   ELSEIF (key .gt. 12 .and. key .le. 15) THEN
      CALL TEP_Primary_Particle(pa, harvest, TEPProd, TEP_n)
   ELSEIF (key .gt. 15 .and. key .le. 17) THEN
      CALL Sedi_Primary_Particle(pa, harvest, sedProd, sedi_n)
   ELSEIF (key .gt. 17 .and. key .le. 20) THEN
      CALL Oil_Primary_Particle(pa, harvest, key, oilProd, oil_n)
   ELSE
      print*, 'in Add_Aggregate:', key, pa%id
   ENDIF
   
   CALL random_number(harvest)   
!   pa%z = 500.0*harvest
   pa%r = pa%pr
   pa%p = pa%Nn * pa%n
  
   i = 1
   age = 1.0d00
   DO WHILE (i .le. 10)
      pa%orgC(i,2) = age*86400.0
      pa%TEP(i,2) = age*86400.0
      age = age * 2.0d00
      i = i + 1
   ENDDO
   i = 1 
   age = 1.0d00
   DO WHILE (i .le. 6)
      pa%oil(i,2) = age*86400.0
      age = age * 2.0d00
      i = i + 1
   ENDDO

   ! Get the size bin of the aggaregate
   i = 1
   CALL find_size(pa%r, i)
   
   ! Define the size bin of the aggregate
   pa%idr = i
   CALL Composite_Stickiness(pa)
   CALL fractal_dimension(pa)
   CALL density(pa)
   CALL velocity(pa)
   !Keep track of the net primary production
   npp(1) = npp(1) + pa%orgC(1,1)*pa%n*mw_orgC                 !g Org_C
   npp(2) = npp(2) + (pa%mineral(1)+pa%mineral(2))*pa%n*mw_co  !g CaCO3
   npp(3) = npp(3) + pa%mineral(3)*pa%n*mw_si                  !g Biogeic Silica
   npp(4) = npp(4) + pa%mineral(4)*pa%n*mw_li                  !g Clay
!   npp(5) = npp(5) + pa%TEP(1,1)*pa%n*mw_TEP                  !g TEP
   npp(5) = npp(5) + pa%TEP(1,1)*pa%n*12.0                     !g TEPC
   npp(6) = npp(6) + pa%oil(1,1)*pa%n                          !g Oil
   npp(7) = npp(7) + pa%sedi(1,1)*pa%n                         !g Sand

   !Why do this again? this is the same as npp?
   inventory(1) = inventory(1) + pa%orgC(1,1)*pa%n
   inventory(2) = inventory(2) + (pa%mineral(1)+pa%mineral(2))*pa%n
   inventory(3) = inventory(3) + pa%mineral(3)*pa%n
   inventory(4) = inventory(4) + pa%mineral(4)*pa%n

END SUBROUTINE Add_Aggregate

!====================================================================

SUBROUTINE OrgC_Primary_Particle(pa, harvest, primaryProd, phyto_n)

   USE the_info
   USE  data_sim
   USE data_stations
   
   IMPLICIT NONE
   TYPE(agg_part), INTENT(INOUT) :: pa
   INTEGER, INTENT(INOUT) :: phyto_n
   DOUBLE PRECISION, INTENT(INOUT) :: harvest, primaryProd
   DOUBLE PRECISION :: calc, arag, opal, clay, orgC
   DOUBLE PRECISION :: vorgC, vcalc, vopal, vclay, tvol
   DOUBLE PRECISION :: rad_p, partVol
   INTEGER :: i

  ! Define different primary particles
   IF (harvest .le. 0.05) THEN
   ! Cocolith particle (mol)		
      pa%orgC(1,1) = 7.d-12
      pa%mineral(1) = 5.d-12
   ELSEIF (harvest .gt. 0.05 .and. harvest .le. 0.07) THEN
   !Aragonite forming phytoplankton particle (mol)
      pa%orgC(1,1) = 7.d-12
      pa%mineral(2) = 5.d-12
   ELSEIF (harvest .gt. 0.07 .and. harvest .le. 0.17) THEN
   !Diatom particle (mol)
      pa%orgC(1,1) = 15.d-12 
      pa%mineral(3) = 5.d-12
   ELSEIF (harvest .gt. 0.17 .and. harvest .le. 1) THEN
   !Picoplankton particle (mol)
      pa%orgC(1,1) = 1.d-12
  ENDIF

!IF (harvest .le. stick_param) THEN
   !Diatom particle (mol)
!print*, 'diatom'
!      pa%orgC(1,1) = 15.d-12 
!      pa%mineral(3) = 5.d-12
!   ELSEIF (harvest .gt. stick_param .and. harvest .le. 0.93) THEN
   !Picoplankton particle (mol)
!print*, 'pico'
!      pa%orgC(1,1) = 1.d-12
!   ELSEIF (harvest .gt. 0.93 .and. harvest .le. 1) THEN
!   ! Cocolith particle (mol)		
!print*, 'cocco'
!      pa%orgC(1,1) = 7.d-12
!      pa%mineral(1) = 5.d-12
!ENDIF




! Calculate the primary particle volume 
   i = 1
   orgC = 0.0
   ! Total OrgC in the primary particle
   DO WHILE (i .le. 10) 
     orgC = orgC + pa%orgC(i,1) 
      i = i + 1
   ENDDO
   ! Total calc, opal and clay in the primary particle
   calc = pa%mineral(1) + pa%mineral(2)
   opal = pa%mineral(3)
   clay = pa%mineral(4)
   ! Change moles to volume
   vorgC = orgC * mw_orgC / rho_orgC
   vcalc = calc * mw_co / rho_co
   vopal = opal * mw_si / rho_si
   vclay = clay * mw_li / rho_li
   ! Get the total volume of the solid material in the primary particle
   tvol = vorgC + vcalc + vopal + vclay
   ! Calculate the radius of the particle
   pa%pr = 1e4*((3.0/(4.0*pi))*tvol)**(1.0/3.0)
   pa%Nn = 1

   ! Calculate the number of particles produced here [mgC/ts
   pa%n = primaryProd/(phyto_n*orgC*mw_orgC*1e3)
!   print*, pa%n, orgC, 'phytoplankton'
   CALL random_number(harvest)
!   pa%z = 1000.0*harvest
   pa%z = 500.0*harvest


CALL particle_volume(pa, tvol)
CALL radius(pa)

END SUBROUTINE OrgC_Primary_Particle

!====================================================================

SUBROUTINE TEP_Primary_Particle(pa, harvest, TEPProd, TEP_n)

    USE the_info
    USE  data_sim
    USE data_stations

    IMPLICIT NONE
    TYPE(agg_part), INTENT(INOUT) :: pa
    INTEGER, INTENT(INOUT) :: TEP_n
    DOUBLE PRECISION, INTENT(INOUT) :: harvest
    DOUBLE PRECISION, INTENT(OUT) :: TEPProd
    DOUBLE PRECISION :: vol, rad, mass

    ! TEP_C [mol] in one particle is:
    pa%TEP(1,1) = 1.0e-12
    ! mass [g] of TEP particle
    mass = pa%TEP(1,1) * mw_TEP
    ! volume [cm3] of TEP particle
    vol = mass / rho_TEP
    ! radius [cm] jof TEP particle
    rad = ((3*vol)/(4*pi))**(1.0/3.0)
    pa%pr = rad * 1e4
    pa%Nn = 1
    ! Calculate the number of particles produced here [mgC/ts
    pa%n = TEPProd/(TEP_n*pa%TEP(1,1)*mw_TEP*1e3)
!   print*, pa%n, pa%TEP(1,1), 'TEP - prod'

    CALL particle_volume(pa, vol)
    CALL radius(pa)
!    pa%z = 1000.0*harvest
    pa%z = 500.0*harvest

END SUBROUTINE TEP_Primary_Particle

!====================================================================

SUBROUTINE Oil_Primary_Particle(pa, harvest, key, oilProd, oil_n)

    ! This subroutine defines a primary oil
    ! particle of defined size

    USE the_info
    USE  data_sim
    USE data_stations

    IMPLICIT NONE
    TYPE(agg_part), INTENT(INOUT) :: pa
    DOUBLE PRECISION, INTENT(INOUT) :: harvest, oilProd
    INTEGER, INTENT(IN) :: key, oil_n

    pa%Nn = 1

!    IF (harvest .le. 0.17) THEN
!        pa%oil(1,1) = 5.0e-10
!        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 0.5
!    ELSEIF (harvest .gt. 0.17 .and. harvest .le. 0.33) THEN
!        pa%oil(1,1) = 5.0e-9
!        pa%n = oilProd/(oil_n*pa%oil(1,1))
!    ELSEIF (harvest .gt. 0.33 .and. harvest .le. 0.5) THEN
!        pa%oil(1,1) = 5.0e-8
!        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 1.5
!    ELSEIF (harvest .gt. 0.5 .and. harvest .le. 0.67) THEN
!        pa%oil(1,1) = 5.0e-7
!        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 2
!    ELSEIF (harvest .gt. 0.67 .and. harvest .le. 0.83) THEN
!        pa%oil(1,1) = 5.0e-6
!        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 4
!    ELSEIF (harvest .gt. 0.83 .and. harvest .le. 1) THEN
!        pa%oil(1,1) = 5.0e-5
!        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 8
!    ELSE
!       print*, 'this should not have happened', pa%id, harvest
!    ENDIF
!
    !Oil droplet mass [g]
    IF (key .eq. 18) THEN
        pa%oil(1,1) = 5.0e-10
        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 0.5
    ELSEIF (key .eq. 19) THEN
        pa%oil(1,1) = 5.0e-9
        pa%n = oilProd/(oil_n*pa%oil(1,1))
    ELSEIF (key .eq. 20) THEN 
        pa%oil(1,1) = 5.0e-8
        pa%n = oilProd/(oil_n*pa%oil(1,1)) * 1.5
    ELSE
       print*, 'this should not have happened', pa%id, key
    ENDIF
    ! Calculate the number density of droplets [#/m^3]
    CALL radius(pa)
!print*, pa%oil(1,1), pa%r, key
!    print*, key, pa%n, (pa%r*1e-3)**3.0*4.0/3.0*3.1415 * pa%n
!    pa%n = 6.25e-6*(pa%pr/1e6)**(-2.3) !1/32xDWH seep3
!    pa%n = 1.25e-5*(pa%pr/1e6)**(-2.3) !1/16xDWH seep2
!    pa%n = 2.5e-5*(pa%pr/1e6)**(-2.3) !1/8xDWH seep
!    pa%n = 5e-5*(pa%pr/1e6)**(-2.3)  !1/4xDWH
!    pa%n = 1e-4*(pa%pr/1e6)**(-2.3) !1/2xDWH
!    pa%n = 2e-4*(pa%pr/1e6)**(-2.3) !DWH - 201-300
!    pa%n = 4e-4*(pa%pr/1e6)**(-2.3) !2xDWH - 201-300 corexit
    pa%o = 1
   CALL random_number(harvest)   
    pa%z = (500 * harvest) !+ stick_param
    !pa%z = 40000.0*harvest + 100000

END SUBROUTINE Oil_Primary_Particle

!====================================================================

SUBROUTINE Sedi_Primary_Particle(pa, harvest, sedProd, sedi_n)

    ! This subroutine defines a primary sedi
    ! particle of defined size

    USE the_info
    USE  data_sim
    USE data_stations

    IMPLICIT NONE
    TYPE(agg_part), INTENT(INOUT) :: pa
    DOUBLE PRECISION, INTENT(INOUT) :: sedProd, harvest
    INTEGER, INTENT(INOUT) :: sedi_n

    !mass of sediment particle [g]
    pa%sedi(1,1) = 5.0e-11 !* stick_param
!    pa%sedi(1,1) = 5.0e-9 !* stick_param  THIS IS TOO BIG
    pa%Nn = 1
    pa%d = 1
    pa%n = sedProd/(pa%sedi(1,1)*sedi_n)

!    pa%n = 1e6
    !Calculate the primary particle sedi mass g

!    CALL particle_volume(pa, partVol)
    CALL radius(pa)
    !pa%z = 1000.0*harvest
    pa%z = 500.0*harvest

END SUBROUTINE Sedi_Primary_Particle
