SUBROUTINE density(p)
	
	! This subroutine defines caculate the density
	! of the aggregate and update the aggregate 
	! properties
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Output variables
	!-----------------
	! update the aggregate properties	
	
	USE the_info
	USE data_sim
   	IMPLICIT NONE

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: rho, porosity, tvol
   	DOUBLE PRECISION :: orgC, calc, opal, clay, dens, TEPC
   	DOUBLE PRECISION :: rho_solid, mass_oil, mass_sedi, rho_sw_m
   	INTEGER :: i, j 

   	! Get the aggragate fractal dimension
   	CALL fractal_dimension(p)
   	! Get the particle volume in cm^3
   	CALL particle_volume(p, tvol)
   	! Get the particle radius in um

   	CALL radius(p)
   	! Get the density of the sea water (measured) g/cm3
   	rho_sw_m = Dens_z(p%idl)

	! Solid masses
   	orgC = 0.0
   	TEPC = 0.0
   	mass_oil = 0.0
   	mass_sedi = 0.0
     ! Get the total organic carbon and TEP
        i = 1
        DO WHILE (i .le. 10)
           orgC = orgC + p%orgC(i,1) 
           TEPC = TEPC + p%TEP(i,1) 
           i = i + 1
        ENDDO
   ! Get the total minerals
   calc = p%mineral(1) + p%mineral(2)
   opal = p%mineral(3)
   clay = p%mineral(4)
   
   	! Get the total oil
   	i = 1
   	DO WHILE (i .le. 6) 
   		mass_oil = mass_oil + p%oil(i,1)
    	i = i + 1
   	 ENDDO   
  
	! Get the total sediment
   	i = 1
   	DO WHILE (i .le. 6) 
   		mass_sedi = mass_sedi + p%sedi(i,1)
    	i = i + 1
   	 ENDDO   
   
	! Get the solid density
	rho_solid = (orgC*mw_orgC/OMCfrac + &
               calc*mw_co + &
               opal*mw_si + &
               clay*mw_li + &
               TEPC*mw_TEP + &
			   mass_oil + &
			   mass_sedi) / tvol
	
	! Commented by Anusha 03/2017
	!IF (p%af .eq. 0 ) THEN
    !  porosity = 1.0 - p%Nn*(p%pr/p%r)**3
    !ELSEIF (p%af .eq. 1 ) THEN
    !  porosity = 0.5
    !ENDIF
    ! Added by Anusha
    ! Define/calculate the porosity
    IF (p%o.eq.1 .or. p%d.eq.1) THEN
		porosity = 0.0
   	ELSEIF (p%af .eq. 0 .and. (p%o.ne.1 .or. p%d.ne.1)) THEN
		porosity = 1.0 - p%Nn*(p%pr/p%r)**3
   	ELSEIF (p%af .eq. 1) THEN
		porosity = 0.5
   	ENDIF
   
   	IF (porosity .lt. 0 .or. porosity .gt. 1) THEN 
		print*, 'porosity less than zero!', porosity, p%frac, p%s, p%Nn, &
        TEPC, tvol, p%af, p%o, p%d
		print*, p%pr, p%r
      	porosity = 0.1
   	ENDIF

	! Get the desity of the porous aggragate
	IF (p%Nn .eq. 1) THEN
		dens = rho_solid
      	porosity = 0
   	ELSEIF (p%Nn .gt. 1) THEN
    	dens = rho_solid*(1.0-porosity) + rho_sw_m*porosity
   	ENDIF
   
	! Update the porosity of the aggregate
   	p%por = porosity
	! Update the density of the aggregate
   	p%rho = dens
!print*, p%id, p%rho, p%r
   	IF (p%por .lt. 0 .or. p%por .gt. 1) THEN
    	   IF (p%Nn .gt. 1) THEN
              print*, 'r lt pr', p%Nn, p%r, tvol
              print*, rho_solid
              print*, p
      	   ENDIF
      	   p%r = p%pr
	   print*, 'p%r = p%pr'
   	ENDIF

   	IF (p%r .eq. 0 .or. p%rho .le. 0) THEN
    	print*, 'in density', p%id, p%n, p%Nn
      	print*, p%rho, p%r, porosity
      	print*, p%orgC(1,1), p%orgC(2,1), p%orgC(3,1)
	ENDIF

   	IF (dens .le. 0) THEN
    	print*, 'rho - problem', porosity, dens, tvol
   	ENDIF
 
   	IF (p%Nn .gt. 1 .and. p%r .eq. p%pr) THEN
    	print*, 'NOT GOOD', p%id, p%Nn, p%r, p%pr, p%frac, p%af, p%o, p%d
   	ENDIF

END SUBROUTINE density
!=====================================================================
SUBROUTINE fractal_dimension_tinna(p)
	
	! This subroutine defines fractal dimension
	! of the aggregate 
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Output variables
	!-----------------
	! update the aggregate fractal demension		
	
	
   	USE the_info
	USE data_sim
	 
   	TYPE(agg_part), INTENT(INOUT) :: p
	
   	IF (p%af .eq. 0) THEN
    	p%frac = 2.0 
   	ELSEIF (p%af .gt. 0) THEN
    	p%frac = 2.999 !2.8
   	ENDIF

END SUBROUTINE fractal_dimension_tinna

SUBROUTINE fractal_dimension(p)
	
	! This subroutine defines fractal dimension
	! of the aggregate 
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Output variables
	!-----------------
	! update the aggregate fractal demension	
	
   	USE the_info
    USE data_sim
	
   	IMPLICIT NONE
   	TYPE(agg_part), INTENT(INOUT) :: p   

	If(p%o .eq. 1 .or. p%d .eq. 1) THEN
		p%frac = 3.0
  	ELSEIF (p%af .eq. 0 .and. p%o .eq. 0 .and. p%d .eq. 0) THEN
    	p%frac = 2.2
   	ELSEIF (p%af .gt. 0) THEN
    	p%frac = 2.999
   	ENDIF
   
END SUBROUTINE fractal_dimension

!=====================================================================
SUBROUTINE dissolve(p)
	
	! This subroutine delete the small aggregates
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	USE the_info
    USE data_sim
   	IMPLICIT NONE

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: orgC, calc, arag, opal, clay, dens, TEPC
   	DOUBLE PRECISION :: mass_oil, mass_sedi, mas_oil_sedi
   	INTEGER :: i, j

   	! Check if the oil and sedi is in the particle
   	mass_oil = 0.0
   	mass_sedi = 0.0
   
   	i = 1
   	DO WHILE (i .le. 6)   
		mass_oil = mass_oil + p%oil(i,1)
	   	i = i + 1
   	ENDDO
   
   	i = 1
   	DO WHILE (i .le. 6)   
		mass_sedi = mass_sedi + p%sedi(i,1)
	   	i = i + 1
   	ENDDO
   
   	mas_oil_sedi = mass_oil + mass_sedi

   	! Find the depth layer of the particle
   	!CALL find_depth(p%z,j)
   	j = p%idl
   
   	IF (p%r .lt. 0.2 .and. p%b .eq. 0 .and. mas_oil_sedi .le. 0.0) THEN
    	i = 1
      	orgC = 0
      	TEPC = 0
      	DO WHILE (i .le. 10) 
        	orgC = orgC + p%orgC(i,1) 
        	TEPC = TEPC + p%TEP(i,1) 
         	i = i + 1
      	ENDDO
      	calc = p%mineral(1) 
      	arag = p%mineral(2)
      	opal = p%mineral(3)
      	clay = p%mineral(4)

      	p%orgC = 0
      	p%mineral = 0
      	p%TEP = 0
  
      	p%b = 1
      	inventory(1) = inventory(1) - orgC*p%n - TEPC*p%n
      	inventory(2) = inventory(2) - (arag+calc)*p%n
      	inventory(3) = inventory(3) - opal*p%n
      	inventory(4) = inventory(4) - clay*p%n

      	lost(1) = lost(1) + orgC*p%n 
      	lost(2) = lost(2) + (arag+calc)*p%n
      	lost(3) = lost(3) + opal*p%n
      	lost(4) = lost(4) + clay*p%n
      	lost(5) = lost(5) + TEPC*p%n
	  
	  	! Keep track of components lost from the aggregate for mass balance 
      	lost_mass_balnce(1) = lost_mass_balnce(1) + orgC*p%n 
      	lost_mass_balnce(2) = lost_mass_balnce(2) + arag*p%n
      	lost_mass_balnce(3) = lost_mass_balnce(3) + calc*p%n	  
      	lost_mass_balnce(4) = lost_mass_balnce(4) + opal*p%n
      	lost_mass_balnce(5) = lost_mass_balnce(5) + clay*p%n
      	lost_mass_balnce(6) = lost_mass_balnce(6) + TEPC*p%n		  

      	orgC_b(j) = orgC_b(j) + orgC*p%n + TEPC*p%n
      	calcC_t(j) = calcC_t(j) + calc*p%n
      	calcA_t(j) = calcA_t(j) + arag*p%n 
      	opal_t(j) = opal_t(j) + opal*p%n 

   	ENDIF

   	i = 1
   	DO WHILE (i .le. 10)
      	!IF (p%orgC(i,1) .lt. 1.0e-50) THEN	! Commented by Anusha 03/2017
		IF (p%orgC(i,1)*p%n .lt. 1.0d-50) THEN  
        	p%orgC(i,1) = 0.0
      	ENDIF
      	i = i + 1
   	 ENDDO

END SUBROUTINE dissolve
!====================================================================

SUBROUTINE disintegrate(p, harvest)

! This subroutine decide the aggragte break up based on
! 1) the particle falls apart if it is not sticky enough to be 
! bound together 
! 2) or, if there is too much shear acting on that particle

! Input variables
!-----------------
! p - pointer for aggregate
! harvest - random number

    USE the_info
    USE data_sim
    IMPLICIT NONE

    TYPE(agg_part), INTENT(INOUT) :: p
    DOUBLE PRECISION, INTENT(INOUT) :: harvest
    DOUBLE PRECISION :: new_Nn, old_Nn, eta

    IF (p%s .lt. 0.02 .and. p%Nn .gt. 1) THEN
       CALL random_number(harvest)
       IF (harvest .gt. 0.99) THEN
          CALL break_p(p, harvest)
       ENDIF
    ENDIF
   
    IF (p%r .gt. 2.0d4 .and. p%Nn .gt. 1) THEN
       CALL random_number(harvest)
       !IF (timestep/86400.0*p%r/1.0e6 .gt. harvest) THEN !Commented by Anusha
       ! get the Kolmogorov scale in the depth layer in cm
       eta = Eta_z(p%idl)
       ! If the particle diameter is greater than eta break up the particle
       IF (p%r*2.0 .gt. eta*10000.0) THEN  
!IF (p%z .gt. 10000 .and. p%r .gt. 10000) print*, 'breaking..r[mm], z[m]', p%r/100, p%z/100, p%id
          CALL break_p(p,harvest)
!IF (p%z .gt. 10000 .and. p%r .gt. 1000) print*, 'broken r[mm], z[m]', p%r/100, p%s
       ENDIF
    ENDIF

END SUBROUTINE disintegrate

!====================================================================
SUBROUTINE particle_volume(p, tvol)
	
	!This subroutine calculte the volume of an aggregate in cm^3
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	! tvol - total volume

	USE the_info
	USE data_sim
   	IMPLICIT NONE
   
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, INTENT(INOUT) :: tvol
   	DOUBLE PRECISION :: ar, avol, orgC, TEPC, calc, opal, clay, mass_oil, mass_sedi
   	DOUBLE PRECISION :: vorgC, vTEPC, vcalc, vopal, vclay, pvol, voil, vsedi
   	INTEGER :: i

   	orgC = 0.0
   	TEPC = 0.0
   	mass_oil = 0.0
   	mass_sedi = 0.0

	! Get the total moles of OrgC and TEP
	i = 1
   	DO WHILE (i .le. 10) 
		orgC = orgC + p%orgC(i,1)
      	TEPC = TEPC + p%TEP(i,1)
      	i = i + 1
   	ENDDO

	! Get the total moles of calcium carbonate
	calc = p%mineral(1) + p%mineral(2)
	! Get total the moles of opal
   	opal = p%mineral(3)
	! Get the total moles of clay
   	clay = p%mineral(4)

	! Get the total mass of oil
   	i = 1
  	DO WHILE (i .le. 6)   
		mass_oil = mass_oil + p%oil(i,1)
	   	i = i + 1
   	 ENDDO

	! Get the total mass of sedi/sediment
   	i = 1
   	DO WHILE (i .le. 6)   
		mass_sedi = mass_sedi + p%sedi(i,1)
	   	i = i + 1
   	ENDDO
  
	! Change moles and mass to volume
   	vorgC = (orgC * mw_orgC / OMCfrac) / rho_orgC
   	vTEPC = TEPC * mw_TEP / rho_TEP
   	vcalc = calc * mw_co / rho_co
  	vopal = opal * mw_si / rho_si
  	vclay = clay * mw_li / rho_li
   	voil = mass_oil / rho_oil   
   	vsedi = mass_sedi / rho_sedi
	 
	! Total volume of the aggregate in cm^3
   	tvol = vorgC + vTEPC + vcalc + vopal + vclay + voil + vsedi

END SUBROUTINE particle_volume

!====================================================================
SUBROUTINE particle_mass(p, mass)
	
	!This subroutine calculte the mass of an aggregate in g
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	! mass - mass of an aggregate
	
  	USE the_info
    USE data_sim
   	IMPLICIT NONE
	
	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, INTENT(INOUT) :: mass
   	DOUBLE PRECISION :: orgC, TEPC, calc, opal, clay 
   	DOUBLE PRECISION :: mass_orgC, mass_oil, mass_sedi
   	INTEGER :: i

   	orgC = 0.0
   	TEPC = 0.0
   	mass_oil = 0.0
   	mass_sedi = 0.0
   
   	! Get the moles of Organic Carbon and TEP
        i = 1
        DO WHILE (i .le. 10) 
           orgC = orgC + p%orgC(i,1) 
           TEPC = TEPC + p%TEP(i,1) 
           i = i + 1
        ENDDO
        mass_orgC = orgC * mw_orgC / OMCfrac
   
   	! Get the total moles of calcium carbonate
   	calc = p%mineral(1) + p%mineral(2)
   	! Get total the moles of opal
   	opal = p%mineral(3)
   	! Get the total moles of clay
   	clay = p%mineral(4)
   
   	! Get the mass of oil (g)
   	i = 1
	DO WHILE (i .le. 6) 
   		mass_oil = mass_oil + p%oil(i,1) 
    	i = i + 1
	ENDDO
   
   	! Get the mass of sediments (g)
    i = 1
 	DO WHILE (i .le. 6) 
    	mass_sedi = mass_sedi + p%sedi(i,1) 
     	i = i + 1
 	ENDDO
	
	! Get the total mass (g)
   	mass = mass_orgC + TEPC * mw_TEP + calc * mw_co + &
          opal * mw_si + clay * mw_li + mass_oil + mass_sedi
   
END SUBROUTINE particle_mass
!====================================================================
SUBROUTINE Composite_Stickiness(p)
	
	! This subroutine defines and calculates the 
	! stickiness of the aggerate based on the stickiness of
	! indvidual components and their volumes
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
   	USE the_info
    USE data_sim
   	IMPLICIT NONE
  
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: stickVar
   	DOUBLE PRECISION :: tvol, TEPC, vTEPC, orgC, vorgC, vbSi
   	DOUBLE PRECISION :: mass_oil, voil, mass_sedi, mass_bSi, vsedi
   	INTEGER :: i

   	CALL particle_volume(p, tvol)

   	orgC = 0.0
   	TEPC = 0.0
   	mass_oil = 0.0
   	mass_sedi = 0.0
   	mass_bsi = 0.0
	
   	i = 1
   	DO WHILE (i .le. 10) 
    	orgC = orgC + p%orgC(i,1)
	TEPC = TEPC + p%TEP(i,1) 
      	i = i + 1
   	ENDDO
 

   	! Get the total mass of oil
   	i = 1
   	DO WHILE (i .le. 6) 
		mass_oil = mass_oil + p%oil(i,1) 
       	i = i + 1
   	ENDDO		

	! Get the total mass of sediment
   	i = 1
   	DO WHILE (i .le. 6) 
		mass_sedi = mass_sedi + p%sedi(i,1) 
       	i = i + 1
   	 ENDDO	   

        mass_bSi = p%mineral(3)
 
   	! Get the volumes of material
   	vTEPC = TEPC * mw_TEP / rho_TEP
   	vorgC = (orgC * mw_orgC / OMCfrac) / rho_orgC
   	voil = mass_oil / rho_oil
   	vsedi = mass_sedi / rho_sedi
   	vbSi = mass_bSi * mw_si / rho_si

	! Different methods to estimate the composite stickiness
	! Pick the one
	!Mtd 1
!	p%s = vTEPC / tvol	
	!Mtd 2
	!p%s = (vTEPC+0.1*vorgC) / tvol	
	!Mtd 3
	!p%s = (vTEPC+1.0*vorgC) / tvol	
	!Mtd 4	
 	p%s =  (vTEPC*s_TEP+ vorgC*s_OrgC + voil*s_oil + vsedi*s_sedi + vbSi*s_bSi) / tvol
!	p%s =  (vTEPC*s_TEP+ vorgC*stick_param + voil*s_oil + vsedi*s_sedi) / tvol
!  	p%s =  (vTEPC*s_TEP+ vorgC*s_OrgC + voil*s_oil + vsedi*stick_param) / tvol
!  	p%s =  (vTEPC*s_TEP+ vorgC*s_OrgC + voil*stick_param + vsedi*s_sedi) / tvol
!  	p%s =  (vTEPC*stick_param + vorgC*s_OrgC + voil*s_oil + vsedi*s_sedi) / tvol
	!Mtd 5
	!p%s =  (vTEPC + vorgC + voil + vsedi) / tvol
	!Mtd 6
	!p%s =  (vTEPC + 0.1*vorgC + voil + vsedi) / tvol	
	
   	IF (p%b .eq. 0) THEN
    	IF (p%s .lt. 0 .or. TEPC .lt. 0) THEN
        	print*, 'stickiness II', p
      	ENDIF
   	ENDIF

END SUBROUTINE Composite_Stickiness
!===================================================================

SUBROUTINE velocity_tinna(p)
	
   	USE the_info
	USE data_sim
   	IMPLICIT NONE
	
	! This subroutine caculate the
	! aggragte terminal velocity in cm/s using the Stokes
	! formulation with White' approximation 
	! details in Jokulsdottir (2016)
	! For oil only aggregates use the functions in tamoc modules
	! dbm_phys
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: d_rho, one, u, fu, du
   	DOUBLE PRECISION :: new_u, d, nu, pm_sign, rho_sw_m
   	INTEGER :: i
		
	! Gett the viscosity of seawater at the aggregate location
	nu = Visc_z(p%idl)
!        nu = 0.01
        p%visc = nu
	! Get the density of the seawater at the aggregate location g/cm3
	rho_sw_m = Dens_z(p%idl)
!        rho_sw_m = 1.023
        p%dens = rho_sw_m
	
   	one = 1
   	i = 1
   	d_rho = ABS(p%rho - rho_sw_m)
   
        CALL find_size(p%r, i)
        p%idr = i
   
   ! Get the sign of the vevlocity
   IF (p%rho .gt. rho_sw_m) THEN
      pm_sign = 1
   ELSEIF (p%rho .lt. rho_sw_m) THEN
      pm_sign = -1
   ENDIF
   
   ! Get the diameter and convert to cm
   d = 2.0*p%r/10000.0
   
   ! Calculate the velocity based on the the Stokes law
   u = g*d_rho*d**2/(18*rho_sw_m*nu)
   
   i = 1
   DO WHILE (one .gt. 0.001 .and. i .lt. 20)
!   DO WHILE (i .lt. 20)
      i = i + 1
      fu = 24*nu*u/d + 6*u**1.5/(1+(d/nu)**0.5) + 0.4*u**2 &
      - 4*g*d_rho*d/(3*rho_sw_m)
      du = 24*nu/d + 9*u**0.5/(1+(d/nu)**0.5) + 0.8*u
      new_u = u - fu/du
      IF (new_u .eq. u) THEN 
         print*, 'problem in velocity', d, p%id
         print*, p
      ENDIF
      one = ABS(new_u - u)
      u = new_u
   END DO
   
   ! Update the aggragate velocity in cm/s
   p%w = pm_sign * u
   IF (p%w .gt. 1e4) print*, 'velocity troubles', p%w, d_rho, d
   IF (p%w .lt. -1e4) print*, 'velocity oops', p%w, d_rho, d

END SUBROUTINE velocity_tinna
!====================================================================
SUBROUTINE velocity_tinna_fix(p)
	
   	USE the_info
	USE data_sim
   	IMPLICIT NONE	
	
	! This subroutine caculate the
	! aggragte terminal velocity in cm/s using the Stokes
	! formulation with White' approximation 
	! details in Jokulsdottir (2016)
	! For oil only aggregates use the functions in tamoc modules
	! dbm_phys
	! Anusha : I got slightly different terms than the ones that
	! Tinna used while implementing Eulers method 
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	
	! Calculates the velocty in cm/s
	! For oil only aggregate and if we want to add gas at some point
	! have functions in tamoc that can be directly taken

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: d_rho, one, u, fu, du
   	DOUBLE PRECISION :: new_u, d, nu, pm_sign, rho_sw_m
	DOUBLE PRECISION :: dnomu, ddnomu, numu, dnumu, coef
   	INTEGER :: i
	
	! Gett the viscosity of seawater at the aggregate location
	nu = Visc_z(p%idl)

	! Get the density of the seawater at the aggregate location g/cm3
	rho_sw_m = Dens_z(p%idl)

   	one = 1
   	i = 1
   	d_rho = ABS(p%rho - rho_sw_m)
   
   ! Get the sign of the velocity
   IF (p%rho .gt. rho_sw_m) THEN
      pm_sign = 1
   ELSEIF (p%rho .lt. rho_sw_m) THEN
      pm_sign = -1
   ENDIF
   
   ! Get the diameter and convert to cm
   d = 2.0*p%r/10000.0
   
   ! Calculate the velocity based on the the Stokes law
   u = g*d_rho*d**2/(18*rho_sw_m*nu)
   coef = 4./3.*d_rho/Dens_z(1)*g*d
   
   DO WHILE (one .gt. 0.001 .and. i .lt. 10)
      i = i + 1
      fu = 24.0*u*nu/d + 6*u**2/(1 + (d*u/nu)**0.5) + 0.4*u**2 - coef
      du = 24.0*nu/d + (12.0*u*nu + 12.0*u**1.5*(nu*d)**0.5 - &
	  3.0*u**1.5*(nu*d)**0.5 + 0.8*u)/((nu)**0.5 + (d*u)**0.5)**2.0 + 0.8*u
	  
      new_u = u - fu/du
      IF (new_u .eq. u) THEN 
         print*, 'problem in velocity', d, p%id
         print*, p
      ENDIF
      one = ABS(new_u - u)
      u = new_u
	  write(*,*) 'u_new', u, fu, du, fu/du
   END DO
   
   ! Update the aggragate velocity in cm/s
   p%w = pm_sign * u
   IF (p%w .gt. 1e4) print*, 'velocity troubles', p%w, d_rho, d
   IF (p%w .lt. -1e4) print*, 'velocity oops', p%w, d_rho, d

END SUBROUTINE velocity_tinna_fix
!====================================================================
SUBROUTINE velocity(p)
	
   	USE the_info
	USE data_sim
   	IMPLICIT NONE	
	
	! This subroutine caculate the
	! aggragte terminal velocity in cm/s using the Stokes
	! formulation with White' approximation and can be used in tratified ambient
	! details in Yick et al (2016)
	! For oil only aggregates use the functions in tamoc modules
	! dbm_phys
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Calculates the velocty in cm/s
	! For oil only aggregate and if we want to add gas at some point
	! have functions in tamoc that can be directly taken

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: d_rho, one, u, fu, du, Nf
   	DOUBLE PRECISION :: u1, du1, u2, du2, coef
   	DOUBLE PRECISION :: new_u, d, half_d, nu, pm_sign
   	DOUBLE PRECISION :: rho_sw_m, rho_a_1, rho_a_2
   	INTEGER :: i, j
	
	! Gett the viscosity of seawater at the aggregate location
	nu = Visc_z(p%idl)
!        nu = 0.01
	
	! Get the density of the seawater at the aggregate location g/cm3
	rho_sw_m = Dens_z(p%idl)
!        rho_sw_m = 1.023
	
	! Density difference between aggregate and ambient water
   	d_rho = ABS(p%rho - rho_sw_m)
!        d_rho = 0.0002
	
	! Get the densities to estimate the buoyancy frequency (s-1)
	if (p%idl.eq.ndz) then
		rho_a_1 = Dens_z(p%idl-1)
		rho_a_2 = Dens_z(p%idl)
	else
		rho_a_1 = Dens_z(p%idl)
		rho_a_2 = Dens_z(p%idl+1)
	endif		
	
	! Calculate the buoynace frequency (1/s)
	if ((rho_a_2-rho_a_1) .le. 0.0 .or. d_rho .lt. 1.0d-7) then
		Nf = 0.0	
	else	
		Nf = (g/rho_a_1*(rho_a_2-rho_a_1)/dz)**0.5
	endif
  !print*, Nf, rho_a_1, nu, '***' 
   	!Get the sign of the velocity
   	IF (p%rho .gt. rho_sw_m) THEN
      pm_sign = 1
   	ELSEIF (p%rho .lt. rho_sw_m) THEN
      pm_sign = -1
   	ENDIF
   
   	!Radius in cm
   	half_d = p%r/10000.0

   	! Calculate the velocity based on the the Stokes law
   	u = 2.0*g*d_rho*half_d**2/(9.0*rho_sw_m*nu)
  
	one = 1
        IF (d_rho .lt. 8e-6) THEN
           one = 0.0
           !print*, 'not dense enough', p%id, p%rho, rho_sw_m
           u = 0.001
        ENDIF
	i = 1  
	
   	DO WHILE (one .gt. 0.001 .and. i .lt. 12)
    	i = i + 1
	  	u1 = 12.0*u*nu/half_d + 6.0*u**2.0/(1.0+(2.0*u*half_d/nu)**0.5) + 0.4*u**2.0
	  	du1 = 12.0*nu/half_d + (12.0*u*nu + 16.9707*u**1.5*(nu*half_d)**0.5 - &
	  	4.2426*u**1.5*(nu*half_d)**0.5)/(nu**0.5 + (2.0*u*half_d)**0.5)**2.0 + 0.8*u
	  	u2 = 1.0 + 1.9*(half_d**3.0*Nf**2.0/(nu*u))**0.5 
	  	du2 = -0.95*(half_d/u)**1.5*Nf/nu**0.5
	  
	  	coef = 8.0/3.0*d_rho/Dens_z(1)*g*half_d

      	fu = u1*u2 - coef
      	du = u1*du2 + u2*du1
      	new_u = u - fu/du
     IF (i .gt. 11) print*, 'VELOCITY',  new_u, u, p%id, p%rho, p%z, p%w
	  	IF (new_u .eq. u) THEN 
        	print*, 'problem in velocity', d, p%id
         	print*, p
      	ENDIF
      IF (new_u .lt. 0) THEN
         print*, 'problem child', p%id, p%r, p%rho
         u = 0.001
         new_u = 0.001
      ENDIF
      	one = ABS(new_u - u)
      	u = new_u
   	END DO

        !*************ATH**********************
        ! if d_rho is almost 0, then W is big
!        IF (d_rho .gt. 0 .and. d_rho .lt. 1e-4) THEN
!           p%w = 0.0
!        !   print*, p%id, d_rho, p%r, p%oil(1,1), p%w*864
!        ENDIF
!        IF (p%w .gt. 10.0) print*, p%id, d_rho, p%r, p%w, '-.-'
!        IF (d_rho .gt. 0 .and. d_rho .lt. 1e-4) print*, "hello!!! really small rho", &
!           2.0*g*d_rho*half_d**2/(9.0*rho_sw_m*nu), p%w
   
   	! Update the aggragate velocity in cm/s
   	p%w = pm_sign * u
        IF (p%w .eq. 0) print*, p
   	IF (p%w .gt. 1e4) print*, 'velocity troubles', p%w, d_rho, d
   	IF (p%w .lt. -1e4) print*, 'velocity oops', p%w, d_rho, d
!   	IF (d_rho .lt. 1e-4 .and. d_rho .gt. 0 .and. p%w*864 .gt. 100) &
!          print*, 'velocity rho', p%w*864, d_rho, p%r, 2.0*g*d_rho*half_d**2/(9.0*rho_sw_m*nu), &
!                  fu/du, one

        ! if single diatom cell
        IF (p%Nn .eq. 1 .and. p%mineral(3) .gt. 0) THEN
           p%w = p%w/10.0
        ENDIF

         ! update in situ viscosity
        p%visc = nu
        p%dens = rho_sw_m

END SUBROUTINE velocity

!====================================================================
SUBROUTINE check_depth(pa, pz, depth, i, j, z_max)
   	
	USE the_info
	USE data_sim	
	
	! This subroutine findout if the particle with index i
	! is wihtin the present depth layer, keep the total count of the 
	! particles added to the depth layer during the count loop and the 
	! indices of all the partcles in the present layer is saved in the 
	! agg_z array

   	TYPE(agg_part), INTENT(INOUT) :: pa
   	DOUBLE PRECISION, INTENT(IN) :: depth
   	INTEGER, INTENT(INOUT) :: pz, j, z_max
   	INTEGER, INTENT(IN) :: i
   	DOUBLE PRECISION :: d

   	d = depth + dz/2.0
   	IF (ABS(pa%z - d) .le. dz/2.0) THEN
		! pointer is assinged the value of the particle index i
	  	! This assigns i to the agg_z(j)
      	pz = i 
	  	! Maximum number of particles in the depth layer 
	  	! in the present count
      	z_max = j
	  	! Keep a count of the particles in the layer
      	j = j + 1
   	ENDIF
	
END SUBROUTINE check_depth
!=====================================================================
SUBROUTINE find_depth(z, i)
	
	! z is depth in cm, i is the depth-box
	! This subroutine finds the depth-box a particle 
	! belongs to when its depth is known

	USE the_info
	USE data_sim
   	DOUBLE PRECISION, INTENT(IN) :: z
   	INTEGER, INTENT(INOUT) :: i
   	INTEGER :: aha

   	i = 1
   	aha = 0
   	DO WHILE (aha .eq. 0) 
    	IF (z .ge. (i-1)*dz .and. z .lt. (i)*dz) THEN
        	aha = 1
      	ELSE 
        	i = i + 1
      	ENDIF
   	 ENDDO

END SUBROUTINE find_depth

!=====================================================================
SUBROUTINE find_size(r, i)

! r is radius in um, i is the radius box

   USE the_info
   USE data_sim

   INTEGER, INTENT(INOUT) :: i
   DOUBLE PRECISION, INTENT(IN) :: r
   INTEGER :: aha

   i = 1
   aha = 0
   DO WHILE (i .lt. n_bins-1 .and. aha .eq. 0)
      if (r .gt. radi(n_bins)) then
         print*, radi(n_bins), r, i, n_bins
         write(*,*) 'radius too large', r, i, n_bins
         stop
      endif
      IF (r .lt. radi(i+1)) THEN
         aha = 1
      ELSE
         i = i + 1
      ENDIF
    ENDDO

   IF (i .gt. n_bins) print*, 'find_size:', i, r

END SUBROUTINE find_size
!=====================================================================
SUBROUTINE check_size(p)
! this subroutine checks if idr is correct

   USE data_sim
   USE the_info

   TYPE(agg_part), INTENT(INOUT) :: p  
   INTEGER :: i, aha
   DOUBLE PRECISION :: r

   i = 1
   aha = 0
   r = p%r
   DO WHILE (i .lt. n_bins-1 .and. aha .eq. 0)
      if (r .gt. radi(n_bins)) then
         print*, radi(n_bins), r, i, n_bins
         write(*,*) 'radius too large', r, i, n_bins
         stop
      endif
      IF (r .lt. radi(i+1)) THEN
         aha = 1
      ELSE
         i = i + 1
      ENDIF
    ENDDO

    IF (p%idr .ne. i) print*, 'BAAAAHHHHH', p%idr, p%r/1e3


END SUBROUTINE check_size
!=====================================================================
SUBROUTINE find_velo(w, i)
	
	! r is depth in um, i is the radius box
	USE the_info
 	USE data_sim
   	INTEGER, INTENT(INOUT) :: i 
   	DOUBLE PRECISION, INTENT(IN) :: w
   	INTEGER :: aha

   	i = 1
   	aha = 0
   	DO WHILE (i .lt. n_bins .and. aha .eq. 0)
      IF (w*864 .lt. radi(i)) THEN
         aha = 1
      ELSE
         i = i + 1
      ENDIF
   	ENDDO

  !IF (i .gt. 60) print*, 'find_velo:', i, w

END SUBROUTINE find_velo
!=====================================================================
SUBROUTINE find_dens(p, i)
	
    USE the_info
    USE data_sim
    IMPLICIT NONE
	
	! This subroutine calculates the density of an aggregate
	! r is depth in um, i is the radius box
	! Input variables
	!-----------------
	! p - pointer for aggregate

   	INTEGER, INTENT(INOUT) :: i 
   	TYPE(agg_part), INTENT(INOUT) :: p  
   	DOUBLE PRECISION :: density, orgC, TEPC, calc, opal, clay
   	DOUBLE PRECISION :: rho_solid, tvol, mass_oil, mass_sedi
   	INTEGER :: aha

   	i = 1
   	aha = 0
   	density = 0.8

   	CALL particle_volume(p, tvol)
	! moles of solid stuff
   	i = 1
   	orgC = 0
   	TEPC = 0
   	DO WHILE (i .le. 10) 
      orgC = orgC + p%orgC(i,1) 
      TEPC = TEPC + p%TEP(i,1) 
      i = i + 1
   	ENDDO
   
   	calc = p%mineral(1) + p%mineral(2)
   	opal = p%mineral(3)
   	clay = p%mineral(4)
   
   	i = 1
   	mass_oil = 0.0
   	DO WHILE (i .le. 6) 
	   mass_oil = mass_oil + p%oil(i,1)
	   i = i + 1
 	ENDDO	
   
    i = 1
    mass_sedi = 0.0
    DO WHILE (i .le. 6) 
 	   mass_sedi = mass_sedi + p%sedi(i,1)
 	   i = i + 1
  	ENDDO	
	
   	rho_solid = (orgC*mw_orgC/OMCfrac + &
               calc*mw_co + &
               opal*mw_si + &
               clay*mw_li + &
               TEPC*mw_TEP + &
			   mass_oil + &
			   mass_sedi) / tvol

   	i = 1
   	DO WHILE (i .le. n_bins .and. aha .eq. 0)
      IF (rho_solid .lt. density) THEN
         aha = 1
      ELSE
         density = density + 0.04
         i = i + 1
      ENDIF
   	ENDDO

   	IF (i .gt. n_bins) THEN
      print*, 'find_dens:', p%r, p%rho
      i = n_bins
   	ENDIF

END SUBROUTINE find_dens
!=====================================================================
SUBROUTINE bottom(p, j)

   !Get the parameters at the bottom of the ocean
   USE the_info
   USE  data_sim
   IMPLICIT NONE

   515 FORMAT (12E15.6E2)

   TYPE(agg_part), INTENT(INOUT) :: p
   INTEGER, INTENT(INOUT) :: j
   DOUBLE PRECISION :: organic, inorganic, TEPC, oil, sedi
   INTEGER :: i, l, m, kk, n

   ! Get the particles at the bottom
   IF (p%z .ge. seafloor) THEN
      p%b = 1
      organic = 0.0
      inorganic = 0.0
      TEPC = 0.0
      oil = 0.0
      sedi = 0.0
  
      i = 1
      DO WHILE (i .le. 10)
         organic = organic + p%orgC(i,1)*p%n ![molC]
         TEPC = TEPC + p%TEP(i,1)*p%n
         seaBed(i) = seaBed(i) + p%orgC(i,1)*p%n + p%TEP(i,1)
         i = i + 1
      ENDDO
	  
	  kk = 1
      DO WHILE (kk .le. 6)
         oil = oil + p%oil(kk,1)*p%n
         kk = kk + 1
      ENDDO
	  
	  kk = 1
      DO WHILE (kk .le. 6)
         sedi = sedi + p%sedi(kk,1)*p%n
         kk = kk + 1
      ENDDO
	  
	  ! Keep track of the orgC at the bottom within the time frame
      sf(1) = sf(1) + organic*mw_orgC	!Org_C 
	  ! Keep track of the CaCO3 at the bottom
      sf(2) = sf(2) + (p%mineral(1)+p%mineral(2))*p%n*mw_co	!CaCO3
	  ! Keep track of the BioSi at the bottom
      sf(3) = sf(3) + p%mineral(3)*p%n*mw_si	!Biogeic Silica
	  ! Keep track of the dust or clay( keyolinite)
      sf(4) = sf(4) + p%mineral(4)*p%n*mw_li	!Clay
	  ! Keep track of the TEP at the bottom
      !sf(5) = sf(5) + TEPC*mw_TEP	!TEP
      sf(5) = sf(5) + TEPC*12           !TEPC [g]
	  ! Keep track of the oil at the bottom
      sf(6) = sf(6) + oil	 
	  ! Keep track of the sedi at the bottom
      sf(7) = sf(7) + sedi	 	   
	  
	  ! Keep track of the total orgC at the bottom
      sf_all(1) = sf_all(1) + organic 
	  ! Keep track of the total TEP at the bottom
      sf_all(2) = sf_all(2) + TEPC	  
	  ! Keep track of the total Mineral 1 at the bottom
      sf_all(3) = sf_all(3) + p%mineral(1)*p%n
	  ! Keep track of the total Mineral 2 at the bottom
      sf_all(4) = sf_all(4) + p%mineral(2)*p%n	  
	  ! Keep track of the total Mineral 3 at the bottom
      sf_all(5) = sf_all(5) + p%mineral(3)*p%n
	  ! Keep track of the total Mineral 4 at the bottom
      sf_all(6) = sf_all(6) + p%mineral(4)*p%n
	  ! Keep track of the total oil at the bottom
      sf_all(7) = sf_all(7) + oil	 
	  ! Keep track of the total sedi at the bottom
      sf_all(8) = sf_all(8) + sedi	  
	  
	  ! Keep track of the different minerals
      seaBed(11) = seaBed(11) + p%mineral(1)*p%n
      seaBed(12) = seaBed(12) + p%mineral(2)*p%n
      seaBed(13) = seaBed(13) + p%mineral(3)*p%n
      seaBed(14) = seaBed(14) + p%mineral(4)*p%n

      inorganic = p%mineral(1)*p%n + p%mineral(2)*p%n

      IF (p%af .gt. 0) THEN
         seaBed(15) = seaBed(15) + organic
         seaBed(16) = seaBed(16) + inorganic
         seaBed(17) = seaBed(17) + p%mineral(3)*p%n
         seaBed(18) = seaBed(18) + p%mineral(4)*p%n
      ENDIF
 
!lets try this instead
      CALL find_size(p%r,j)

      sizeBed(j) = sizeBed(j) + organic   !mol orgC/m2

      CALL find_velo(p%w, l)
      veloBed(l) = veloBed(l) + organic

      m = 0
      CALL find_dens(p, m)
      densBed(m) = densBed(m) + organic
  
     ! Track the oil settling at the bottom
      oilBed(j) = oilBed(j) + oil
  
       ! Track the sediments settling at the bottom
      sediBed(j) = sediBed(j) + sedi  
    
      WRITE(115,515) p%Nn, p%r, p%w, p%rho, organic, inorganic, sedi, p%por, p%pr, oil, p%mineral(3)*p%n, TEPC
  
   ENDIF

END SUBROUTINE bottom
!=====================================================================

SUBROUTINE pool(p1, p2)
	
    USE the_info
    USE data_sim
    IMPLICIT NONE	
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Here two smaller particles are put togehter and saved as one and 
	! the other is made a particle at the bottom

   	TYPE(agg_part), INTENT(INOUT) :: p1, p2
   	INTEGER :: i

   	i = 1
   	DO WHILE (i .le. 10) 
   	 	p1%orgC(i,1) = (p1%n*p1%orgC(i,1) + p2%n*p2%orgC(i,1))/(p1%n+p2%n)
      	p1%TEP(i,1) = (p1%n*p1%TEP(i,1) + p2%n*p2%TEP(i,1))/(p1%n+p2%n)
      	i = i + 1
   	ENDDO
             
   	i = 1
   	DO WHILE (i .le. 4) 
    	p1%mineral(i) = (p1%n*p1%mineral(i) + p2%n*p2%mineral(i))/(p1%n+p2%n)
      	i = i + 1
   	ENDDO
   
   	i = 1
   	DO WHILE (i .le. 6)
		p1%oil(i,1) = (p1%n*p1%oil(i,1) + p2%n*p2%oil(i,1))/(p1%n+p2%n)
       	i = i + 1
    ENDDO
	
    i = 1
    DO WHILE (i .le. 6)
 		p1%sedi(i,1) = (p1%n*p1%sedi(i,1) + p2%n*p2%sedi(i,1))/(p1%n+p2%n)
        i = i + 1
    ENDDO	
	
   	p1%n = p1%n + p2%n
   	p1%p = p1%p + p2%p
   	p1%z = (p1%z*p1%n + p2%z*p2%n)/(p1%n+p2%n)
   
   	! Make the particle p2 a bottom particle
   	p2%b = 1

END SUBROUTINE pool
!=====================================================================

SUBROUTINE radius_tinna(p)
	
    USE the_info
    USE data_sim
    IMPLICIT NONE
	
	! This subroutine caclculate the aggreagte
	! radius
   	! this subroutine is called by stick() and break_p and respiration
   	! calculates new radius after coagulation or fragmentation
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: tvol, pvol, vol_fp

   	CALL fractal_dimension(p)
	CALL particle_volume(p, tvol)

   	pvol = tvol/p%Nn

   	p%pr = ((3.0*pvol/(4.0*pi))**(1.0/3.0)*1e4)
   	if (p%pr .gt. 1000.) then
   		write(*,*) 'p%pr',p%pr, tvol, pvol, p%Nn
	endif
   	IF (p%af .eq. 0) THEN !if particle is an aggregate:
    	p%r = p%pr*p%Nn**(1.0/p%frac)
   ELSEIF(p%af .eq. 1) THEN !particle is fecal pellet:
      vol_fp = tvol * 2.0
      p%r = ((3.0*vol_fp/(4.0*pi))**(1.0/3.0)*1e4)
   ELSE 
      p%r = p%pr*p%Nn**(1.0/p%frac)
   ENDIF

END SUBROUTINE radius_tinna

SUBROUTINE radius(p)
   	USE the_info
	USE data_sim
	IMPLICIT NONE
	! This subroutine caclculate the aggreagte
	! radius
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Get the particle radius in um
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: tvol, pvol, vol_fp, r0
	INTEGER :: i

	! Save the intial aggragate size
	r0 = p%pr
      

	! Get the aggragate fractal radius
   	CALL fractal_dimension(p)
	! Get the particle volume in cm3
   	CALL particle_volume(p, tvol)
	
   	! Calculate the primary particle volume  
   	pvol = tvol/p%Nn


   	! Calculate the primary particle radius um
   	p%pr = 1.0e4*((3.0/(4.0*pi))*pvol)**(1.0/3.0)

   	! Check if the particle is all oil then it is not an aggregate
	IF (p%o .eq. 1 .or. p%d .eq. 1) THEN 
		! Particle aggragate is same as the compact radius
		p%r = ((3.0*tvol/(4.0*pi))**(1.0/3.0)*1e4)
		
	! Check if the particle is an aggragate	
   	ELSEIF (p%af .eq. 0 .and. p%o .eq. 0 .and. p%d .eq. 0) THEN 
		! Then particle is an aggregate
      	p%r = p%pr*p%Nn**(1.0/p%frac)

   	! Particle is fecal pellet:
   	ELSEIF(p%af .gt. 0) THEN
   	 	!why multiplied by 2 - May be because porosity is defined 
		!as 0.5 for fecal pellets
		vol_fp = tvol * 2.0		
      	p%r = ((3.0*vol_fp/(4.0*pi))**(1.0/3.0)*1e4)
	  
   	ENDIF

   ! Get the radius size bin of the aggregate 
   if (p%r.ge.r0) then
	   i = p%idr
   else
	   i = 1
   endif
  IF (p%r .gt. 1e6) print*, 'radius is pretty big', p%id, p%r, p%s
   CALL find_size(p%r, i)
   ! Update the radius size bin of the aggregate
   p%idr = i

END SUBROUTINE radius

!=====================================================================
SUBROUTINE collide4(p1, p2, harvest)
   	USE the_info
	USE data_sim
   	IMPLICIT NONE
	
	! This subroutine checks the collision
	! of different aggreates
	
	! Input variables
	!-----------------
	! p1, p2 - pointer for aggregate

   	TYPE(agg_part), INTENT(INOUT) :: p1, p2
   	DOUBLE PRECISION, INTENT(INOUT) :: harvest
   	DOUBLE PRECISION :: test1, test2, endofdt, semiTime
   	DOUBLE PRECISION :: q, dt, p, muu, nu, epsi, rho_sw_m
   	DOUBLE PRECISION :: beta, betaBr, betaSh, betaDS
   	DOUBLE PRECISION :: r1, r2, w1, w2, prob, alpha
   	DOUBLE PRECISION :: organic, inorganic, dwTotal, s, dw_orgC
   	DOUBLE PRECISION :: r, mean, meanp
   	INTEGER :: q_int

   	r1 = p1%r/1e6  ! r1 in m
   	r2 = p2%r/1e6  ! r2 in m
   	w1 = p1%w/100  ! w1 in m/s
   	w2 = p2%w/100  ! w1 in m/s
 
   	p = MIN(r1,r2)/MAX(r1,r2)
   
   	!turbulent energy dissipation rate [cm2/s3]
	epsi = Epsilon_z(p1%idl)
		
	!turbulent energy dissipation rate [m2/s3] 
	!epsi = 1e-8
	epsi = epsi/1.0e4
	
   	! Get the dynamic viscosity at the depth level g/cms
   	muu = Visc_z(p1%idl)
	
	!Get the density of the sea water at the depth level g/m^3
	rho_sw_m = Dens_z(p1%idl)
   
   	!Calculate nu kinematic viscosity [m2/s]
   	nu = (muu/10.0)/(rho_sw_m*1000.0)                  

   	!Get the total number of primary particles before coagulation 
	test1 = p1%p + p2%p

	!Set the time limit for the loop
   	endofdt = timestep
   
   	! Make the time step the maximum possible (ie equal to the time limit)
   	dt = timestep
   
   	! Variable to keep track of the current time witin the loop 
   	semiTime = 0
   
   	! Total number of primary particles in the two colliding particles 
	! again?? Why?
   	test1 = p1%p + p2%p
   
   	! Calculate the rate of coagulation and the porbability of aggregation 
   	DO WHILE (semiTime .lt. endofdt .and. MAX(p1%Nn,p2%Nn) .lt. 1e9)
		
	   	r1 = p1%r/1e6  ! r1 in m
	   	r2 = p2%r/1e6  ! r2 in m
	   	w1 = p1%w/100  ! w1 in m/s
	   	w2 = p2%w/100  ! w1 in m/s
 
	   	p = MIN(r1,r2)/MAX(r1,r2)		
	   
	   ! Get the stickiness of the particles
      alpha = MAX(p1%s,p2%s)
	  
	  !Brownian motion:
	  ! Commented by Anusha 05/2017
      !betaBr = 8*k*T/(6*mu)*(r1+r2)**2/(r1*r2)   ! in m^3/s
	  betaBr = 8*k*(Temp_z(p1%idl)+273.15)/(6*Visc_z(p1%idl)/10.0)*(r1+r2)**2/(r1*r2)
!	  betaBr = 8*k*(Temp_z(p1%idl)*stick_param+273.15)/(6*Visc_z(p1%idl)/10.0)*(r1+r2)**2/(r1*r2)
	  
	  !Shear (Hill, 1992):
     !IF (p1%z .le. 10000) THEN
        !betaSh = 4/3.0*gamma*(r1+r2)**3              ! in m^3/s Rectil.
     !ELSEIF (p1%z .gt. 10000) THEN
        !betaSh = 4/3.0*(gamma/100)*(r1+r2)**3          ! in m^3/s
     !ENDIF
      betaSh = 9.8*(p**2/(1+2*p)**2)*(epsi/nu)**0.5*(r1+r2)**3  ! in m^3/s C
     !betaSh = 9.8*(p**2/(1+2*p**2))*(epsi/nu)**0.5*(r1+r2)**3  ! in m^3/s C
	  
	  !Differential settling
      !betaDS = pi*(r1+r2)**2*ABS(w2-w1)          ! rectilinear
      betaDS = 0.5*pi*MIN(r1,r2)**2*ABS(w2-w1)    ! in m^3/s

	  ! Get the coagulation kernal as the addition of Brownian motion,
	  ! Shear and Differential settling
      beta = betaBr + betaSh + betaDS                   ! in m^3/s
	  
	  !Probability of aggregation
      prob = alpha*beta*dt*100.0/dz   !*MIN(p1%n,p2%n)	! why there is no p1%n and p2%n??
	  ! Use the modified Euler??
	  ! Check to see if the time probability time step is larger
	  ! this is decided yes of prob>1. If so reduce the time step used 
	  ! to calculate the prob. So the reduction of dt is done in the following
	  ! loop  as explained in Jokulsdottir (2016).
      DO WHILE (prob .gt. 1.0) 
!print*, '***************************', prob
         dt = dt/10.0
         prob = alpha*beta*dt*100.0/dz!*MIN(p1%n,p2%n)
      ENDDO

	  ! Find number of j-particles to collide with i-particles
	  ! Equation 17 of Jokulsodottir
      q = MAX(p1%n,p2%n)/MIN(p1%n,p2%n)
      mean = q * (1 - (1-prob)**MIN(p1%n,p2%n))
      meanp = mean - FLOOR(mean)
!	 print*, mean, q, prob, dt
      CALL random_number(harvest)
      IF (meanp .gt. harvest) THEN
         q_int = CEILING(mean)
      ELSEIF(meanp .lt. harvest) THEN
         q_int = FLOOR(mean)
      ENDIF

      IF (q_int .lt. 0) THEN
         print*, 'ERROR: collide: q_int lt 0', q, mean, q_int, p1%id, p2%id, p1%n, p2%n
         q_int = 0
      ENDIF
	 
      IF (q_int .gt. 0) THEN
         CALL stick2(p1, p2, q_int)
      ENDIF

      semiTime = semiTime + dt
   ENDDO

	! Keep track of the mass of the total number of primary particles after coagulation
   test2 = p1%p + p2%p

   IF (ABS(test1-test2) .gt. 10) THEN
      print*, 'ERROR: collide4', test1, test2, p1%id, p2%id, p1%Nn, p2%Nn, p1%n, p2%n
   ENDIF

!  IF (p1%id .eq. 24 .or. p2%id .eq. 24) THEN
!  IF (p1%r .gt. 8e5 .or. p2%r .gt. 8e5) THEN
!     print*, 'radius is big', p1%id, p2%id
!     print*, p1%s, p2%s, p1%Nn, p2%Nn
!     print*, p1%n, p2%n, p1%r, p2%r
!     print*, p1%orgC(1,1), p2%orgC(1,1), p1%TEP(1,1), p2%TEP(1,1)
!     print*, p1%sedi(1,1), p2%sedi(1,1), p1%oil(1,1), p2%oil(1,1)
!     print*, p1%w, p2%w, p1%z, p2%z
!     print*, '-----------------big---------------------'
!  ENDIF
 
  IF (p1%oil(1,1) .lt. 0.0 .or. p2%oil(1,1) .lt. 0.0) THEN
     print*, p1%id, p2%id, p1%oil(1,1), p2%oil(1,1)
  ENDIF

END SUBROUTINE collide4
!========================================================================
SUBROUTINE stick2(p1, p2, intRate)
	
    USE the_info
    USE data_sim	
	
	! This subroutine estimated the amount of material in the aggregates
	! after the collision of different types of aggreates
	
	! Input variables
	!-----------------
	! p - pointer for aggregate

   	TYPE(agg_part), INTENT(INOUT) :: p1, p2
   	INTEGER, INTENT(INOUT) :: intRate
   	INTEGER :: i, whereAmI, if1, if2, if3
   	DOUBLE PRECISION :: test1, test2, x, bak_n1, bak_n2, z_mid
   	DOUBLE PRECISION :: tot_TEP, tot_OrgC, tot_miral, tot_sedi, tot_oil
   	DOUBLE PRECISION :: None_oil_mass, None_sedi_mass
   
   	! Flags to keep track of collision events
   	if1 = 0
   	if2 = 0
   	if3 = 0

   	! Keep track of the mass of the total number of primary particles before coagulation
   	test1 = p1%orgC(1,1)*p1%n + p2%orgC(1,1)*p2%n
   	! what happens if only TEP or only oil particles collide they do not have OrgC then 
   	! this check if not done?    
   
   	x = 2.0  !this is to get out of problems with TYPE and MOD

   	bak_n1 = p1%n
   	bak_n2 = p2%n
   
   	! Check for number of collisions where p%n can go below 1
   	IF (p1%n .eq. p2%n .and. p1%n .le. 1) THEN
		n_neglect = n_neglect + 1
   	ENDIF
   
   	IF (p1%n .gt. p2%n .and. intRate*p2%n .lt. p1%n) THEN
		if1 = 1
      	whereAmI = 1
	  	! Get the new number of primary pariticles in one aggregate in p2
      	p2%Nn = p1%Nn*intRate + p2%Nn

	  	! Initialize ALD
	  	tot_TEP = 0.0
	  	tot_OrgC = 0.0
	  	tot_miral = 0.0
	  	tot_sedi = 0.0
	  	tot_oil = 0.0
	 
      	i = 1
      	DO WHILE (i .le. 10)
        	p2%orgC(i,1) = p2%orgC(i,1) + intRate*p1%orgC(i,1)
         	p2%TEP(i,1) = p2%TEP(i,1) + intRate*p1%TEP(i,1)
		 	tot_TEP = tot_TEP + p2%TEP(i,1)
		 	tot_OrgC = tot_OrgC + p2%orgC(i,1)
         	i = i + 1
      	 ENDDO
	  
      	i = 1
      	DO WHILE (i .le. 4)
        	p2%mineral(i) = p2%mineral(i) + intRate*p1%mineral(i)
		  	tot_miral = tot_miral + p2%mineral(i)		  
          	i = i + 1
      	ENDDO

	  	! Adjust the oil in the particles
      	i = 1
      	DO WHILE (i .le. 6)
			p2%oil(i,1) = p2%oil(i,1) + intRate*p1%oil(i,1)
		  	tot_oil = tot_oil + p2%oil(i,1)		  
		  	i = i + 1
	  	ENDDO
	  
	  	! Adjust the sedi in the particles
      	i = 1
      	DO WHILE (i .le. 6)
			p2%sedi(i,1) = p2%sedi(i,1) + intRate*p1%sedi(i,1)
		  	tot_sedi = tot_sedi + p2%sedi(i,1)
		  	i = i + 1
	  	ENDDO
	
    	None_oil_mass = tot_TEP + tot_OrgC + tot_miral + tot_sedi
    	! Check if the particle is only oil
   		IF (p2%o .eq. 1) THEN
    		IF (None_oil_mass .ne. 0.0 .or. p2%Nn .gt.1) THEN
 				p2%o = 0
 			ENDIF
 		ENDIF
	
    	None_sedi_mass = tot_TEP + tot_OrgC + tot_miral + tot_oil
   	 	! Check if the particle is only sedi
   		IF (p2%d .eq. 1) THEN
    		IF (None_sedi_mass .ne. 0.0 .or. p2%Nn .gt.1) THEN
 				p2%d = 0
 			ENDIF
 		ENDIF	
	
		! After collision place the new aggregates in the middle of 
		! thier depth levels 
                ! Would be better if this was a weighted average
		z_mid = (p1%z + p2%z)/2.
		p1%z = z_mid
		p2%z = z_mid

      	p1%n = p1%n - intRate*p2%n
	  
	  	!Update the aggregate radius
      	CALL radius(p2)
      	CALL density(p2)
	  	CALL Composite_Stickiness(p2)
	  	CALL velocity(p2)
!	  	CALL velocity_tinna(p2)
		
	ELSEIF (p1%n .lt. p2%n .and. intRate*p1%n .lt. p2%n) THEN
		if2 = 1
      	whereAmI = 2
	  	! Get the new number of primary particles in p1 after aggregation
      	p1%Nn = p1%Nn + p2%Nn*intRate
	  
	  	!Initialize
	  	tot_TEP = 0.0
	  	tot_OrgC = 0.0
	  	tot_miral = 0.0
	  	tot_sedi = 0.0
	  	tot_oil = 0.0
	  
      	i = 1
      	DO WHILE (i .le. 10)
        	p1%orgC(i,1) = p1%orgC(i,1) + intRate*p2%orgC(i,1)
         	p1%TEP(i,1) = p1%TEP(i,1) + intRate*p2%TEP(i,1)
		 	tot_TEP = tot_TEP + p1%TEP(i,1)
		 	tot_OrgC = tot_OrgC + p1%orgC(i,1)
         	i = i + 1
      	ENDDO
	  
      	i = 1
      	DO WHILE (i .le. 4)
        	p1%mineral(i) = p1%mineral(i) + intRate*p2%mineral(i)
		 	tot_miral = tot_miral + p1%mineral(i)		 
         	i = i + 1
      	ENDDO
	  
      	i = 1
      	DO WHILE (i .le. 6)
	  		p1%oil(i,1) = p1%oil(i,1) + intRate*p2%oil(i,1)
			tot_oil = tot_oil + p1%oil(i,1)
        	i = i + 1
     	ENDDO
  
     	i = 1
     	DO WHILE (i .le. 6)
  			p1%sedi(i,1) = p1%sedi(i,1) + intRate*p2%sedi(i,1)
			tot_sedi = tot_sedi + p1%sedi(i,1)
       		i = i + 1
    	ENDDO

    	None_oil_mass = tot_TEP + tot_OrgC + tot_miral + tot_sedi
    	
		! Check if the particle is only oil
   		IF (p1%o .eq. 1) THEN
    		IF (None_oil_mass .ne. 0.0 .or. p1%Nn .gt.1) THEN
 				p1%o = 0
 			ENDIF
 		ENDIF	
	
    	None_sedi_mass = tot_TEP + tot_OrgC + tot_miral + tot_oil
    	! Check if the particle is only oil
   		IF (p1%d .eq. 1) THEN
    		IF (None_sedi_mass .ne. 0.0 .or. p1%Nn .gt.1) THEN
 				p1%d = 0
 			ENDIF
 		ENDIF	
	
		! After collision place the new aggregates in the middle of 
		! thier depth levels 
		z_mid = (p1%z + p2%z)/2.
		p1%z = z_mid
		p2%z = z_mid
	  
	  	! Update the aggregate properties
      	p2%n = p2%n - intRate*p1%n
      	CALL radius(p1)
      	CALL density(p1)
	  	CALL Composite_Stickiness(p1)
	  	CALL velocity(p1)
!	  	CALL velocity_tinna(p1)
	  
   	 ELSEIF (p1%n .eq. p2%n .and. p1%n .gt. 1) THEN
	 	if3 = 1               
      	whereAmI = 3
      	p1%Nn = p1%Nn + p2%Nn
      	p2%Nn = p1%Nn
 
      	i = 1
	  	!Initialize
	  	tot_TEP = 0.0
	  	tot_OrgC = 0.0
	  	tot_miral = 0.0
	  	tot_sedi = 0.0
	  	tot_oil = 0.0
	  
      	DO WHILE (i .le. 10)
        	p1%orgC(i,1) = p1%orgC(i,1) + p2%orgC(i,1) 
         	p2%orgC(i,1) = p1%orgC(i,1)
         	p1%TEP(i,1) = p1%TEP(i,1) + p2%TEP(i,1) 
         	p2%TEP(i,1) = p1%TEP(i,1)
		 	tot_TEP = tot_TEP + p1%TEP(i,1)
		 	tot_OrgC = tot_OrgC + p1%orgC(i,1)
         	i = i + 1
      	ENDDO
	  
      	i = 1
      	DO WHILE (i .le. 4)
        	p1%mineral(i) = p1%mineral(i) + p2%mineral(i)  
         	p2%mineral(i) = p1%mineral(i)
		 	tot_miral = tot_miral + p1%mineral(i)
         	i = i + 1
      	 ENDDO
	  
      	i = 1
      	DO WHILE (i .le. 6)
	  		p1%oil(i,1) = p1%oil(i,1) + p2%oil(i,1) 
	  		p2%oil(i,1) = p1%oil(i,1)
			tot_oil = tot_oil + p1%oil(i,1)
        	i = i + 1
     	ENDDO
  
     	i = 1
     	DO WHILE (i .le. 6)
  			p1%sedi(i,1) = p1%sedi(i,1) + p2%sedi(i,1) 
  			p2%sedi(i,1) = p1%sedi(i,1)
			tot_sedi = tot_sedi + p1%sedi(i,1)
       		i = i + 1
    	ENDDO
		  
      	IF (MOD(p1%n,x) .eq. 0.0) THEN 
        	p1%n = p1%n/2.0
         	p2%n = p1%n
      	ELSE
        	p1%n = (p1%n+1)/2.0
         	p2%n = (p2%n-1)/2.0
      	ENDIF
	  
    	! Check if the particle is only oil
   	 	None_oil_mass = tot_TEP + tot_OrgC + tot_miral + tot_sedi
    	! Check if the particle is only oil
   		IF (p1%o .eq. 1) THEN
    		IF (None_oil_mass.ne. 0.0 .or. p1%Nn .gt.1) THEN
 				p1%o = 0
 			ENDIF
 		ENDIF
	
    	None_sedi_mass = tot_TEP + tot_OrgC + tot_miral + tot_oil
    	! Check if the particle is only oil
   		IF (p1%d .eq. 1) THEN
    		IF (None_sedi_mass.ne. 0.0 .or. p1%Nn .gt.1) THEN
 				p1%d = 0
 			ENDIF
 		ENDIF
	
		! After collision place the new aggregates in the middle of 
		! thier depth levels 
		z_mid = (p1%z + p2%z)/2.
		p1%z = z_mid
		p2%z = z_mid

      	CALL radius(p1)
      	CALL radius(p2)
      	CALL density(p1)
      	CALL density(p2)
	  	CALL Composite_Stickiness(p1)
	  	CALL Composite_Stickiness(p2)
	  	CALL velocity(p1)
	  	CALL velocity(p2)
!	  	CALL velocity_tinna(p1)
!	   	CALL velocity_tinna(p2)
	
	ENDIF
	   
   ! Calcualte the number of primary particles
   p1%p = p1%Nn * p1%n
   p2%p = p2%Nn * p2%n

	! Keep track of the mass of the total number of primary particles after coagulation
   test2 = p1%orgC(1,1)*p1%n + p2%orgC(1,1)*p2%n

   ! Check for the conservation of mass of the total number of primary particles
   IF (ABS(test1-test2) .gt. 1)  print*, 'stick2: failed conservation of mass', test1, test2, rate

   ! Tracking for un expected descrepencies??
   IF (p1%n .lt. 0) print*, 'p1 less than 0', p1%Nn, p2%Nn, WhereAmI, &
      bak_n1, bak_n2, p1%n, p2%n, intRate, p1%id, p2%id
   IF (p2%n .lt. 0) print*, 'p2 less than 0', p1%Nn, p2%Nn, WhereAmI, &
      bak_n1, bak_n2, p1%n, p2%n, intRate, p1%id, p2%id
   IF (p1%orgC(1,1) .lt. 0 .or. p2%orgC(1,1) .lt. 0) print*, 'orgC!!!', &
      WhereAmI, p1%orgC(1,1), p1%Nn, p1%n, p2%Nn, p2%n, intRate
   IF (p1%Nn .lt. 0) print*, 'N1 lt 0', WhereAmI, rate, bak_n1, bak_n2
   IF (p2%Nn .lt. 0) print*, 'N2 lt 0', WhereAmI, rate, bak_n1, bak_n2

	n_collision = n_collision + 1
	if (if1.eq.0 .and. if2.eq.0 .and. if3.eq.0) then
		n_stick = n_stick + 1
	endif

END SUBROUTINE stick2
!========================================================================
SUBROUTINE self_collide(p, harvest)
	
    USE the_info
 	USE data_sim
	
	! This subroutine checks the collision
	! of similar aggreates

	! Input variables
	!-----------------
	! p - pointer for aggregate
	
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, INTENT(INOUT) :: harvest
   	DOUBLE PRECISION :: beta, betaBr, betaSh, betaDS, del_t
   	DOUBLE PRECISION :: r1, w1,prob, alpha
   	DOUBLE PRECISION :: test1, test2, gamma_0

   	! Keep track of the total number of primary particles before coagulation
   	test1 = p%p

   	r1 = p%r/1.0e6  ! r1 in m
   	gamma_0 = Gamma_z(p%idl)

   IF (p%n .gt. 2) THEN
      alpha = p%s
	  betaBr = 8*k*(Temp_z(p%idl)+273.15)/(6*Visc_z(p%idl)/10.0)*(2*r1)**2/(r1**2)
!	  betaBr = 8*k*(Temp_z(p%idl)*stick_param+273.15)/(6*Visc_z(p%idl)/10.0)*(2*r1)**2/(r1**2)
      IF (p%z .le. 10000) THEN
         betaSh = 4/3.0*gamma_0*(r1)**3              ! in m^3/s
      ELSEIF (p%z .gt. 10000) THEN
         betaSh = 4/3.0*(gamma_0/100)*(r1)**3          ! in m^3/s
      ENDIF
      betaDS = 0.0
      beta = (betaBr + betaSh + betaDS)
      !prob = 1 - exp(-alpha*beta*sf1*timestep)
      prob = alpha*beta*p%n*timestep


      CALL random_number(harvest)
      IF (prob .gt. harvest) THEN
         CALL self_stick(p)
      ENDIF
   ENDIF

   IF (prob .gt. 1) THEN
   		print*, 'self_collide', p%id, p%Nn, p%n
   ENDIF
   
   ! Keep track of the total number of primary particles after coagualation
   test2 = p%p
   IF (ABS(test1-test2) .gt. 1) THEN 
      print*, 'self', p%id, p%Nn
   ENDIF

END SUBROUTINE self_collide
!========================================================================
!========================================================================
SUBROUTINE self_stick(p)
	
    USE the_info
  	USE data_sim	
	
	! This subroutine estimated the amount of material in the aggregates
	! after the collision of similar aggreates
	
	! Input variables
	!-----------------
	! p - pointer for aggregate

   	TYPE(agg_part), INTENT(INOUT) :: p
   	INTEGER :: i
   	DOUBLE PRECISION :: test1, test2, tot_TEP, tot_OrgC, tot_miral, tot_sedi, tot_oil
   	DOUBLE PRECISION :: None_oil_mass, None_sedi_mass
   
   	! Initialize parameters
   	! Total TEP in the aggregate
   	tot_TEP = 0.0
   	! Total OrgC in the aggregate
   	tot_OrgC = 0.0
   	! Total mineral in the aggregate
   	tot_miral = 0.0
   	! Total sedi in the aggregate
   	tot_sedi = 0.0
   	! Total oil in the aggregate
   	tot_oil = 0.0   
   
   	! Keep track of the mass of the total primary particles 
   	! before coagulation
   	test1 = p%orgC(1,1)*p%n
 
   	p%n = p%n / 2.0
   	p%Nn = p%Nn * 2.0
   
   	i = 1
   	DO WHILE (i .le. 4)
    	p%mineral(i) = p%mineral(i) * 2.0
	  	tot_miral = tot_miral + p%mineral(i)
      	i = i + 1
   	ENDDO

   	i = 1 
   	DO WHILE (i .le. 10) 
    	p%orgC(i,1) = p%orgC(i,1) * 2.0
      	p%TEP(i,1) = p%TEP(i,1) * 2.0
	  	tot_TEP = tot_TEP + p%TEP(i,1)
	  	tot_OrgC = tot_OrgC + p%orgC(i,1)
      	i = i + 1
   	 ENDDO
   
   	! Get the total oil in the particle
   	i = 1 
   	DO WHILE (i .le. 6)
		p%oil(i,1) = p%oil(i,1)*2.0
	   	tot_oil = tot_oil + p%oil(i,1)
	   	i = i + 1
   	ENDDO		
   
   	! Get the total sedi in the particle
   	i = 1 
   	DO WHILE (i .le. 6)
		p%sedi(i,1) = p%sedi(i,1)*2.0
	   	tot_sedi = tot_sedi + p%sedi(i,1)
	   	i = i + 1
  	ENDDO
     
   	None_oil_mass = tot_TEP + tot_OrgC + tot_miral + tot_sedi
   	! Check if the particle is only oil
  	IF (p%o .eq. 1) THEN
   		IF (None_oil_mass.ne. 0.0) THEN
			p%o = 0
		ENDIF
	ENDIF
	
    None_sedi_mass = tot_TEP + tot_OrgC + tot_miral + tot_oil
    ! Check if the particle is only oil
   	IF (p%d .eq. 1) THEN
    	IF (None_sedi_mass .ne. 0.0) THEN
 			p%d = 0
 		ENDIF
 	ENDIF	

	!Update the aggregate properties
   	CALL radius(p)
        CALL check_size(p)
   	CALL density(p)
   	CALL Composite_Stickiness(p)
   	CALL velocity(p)
!	CALL velocity_tinna(p)
   	! Keep track of the mass of the total number of primary particles
   	! after coagulation
   	test2 = p%orgC(1,1)*p%n

   	IF (ABS(test1-test2) .gt. 1) THEN 
    	print*, 'self_stick', p%Nn, p%n
   	ENDIF
	print*, 'SELF-STICK', p%id, p%TEP(1,1)

END SUBROUTINE self_stick
!========================================================================
SUBROUTINE break_p(p, harvest)
	
    USE the_info
    USE data_sim
    IMPLICIT NONE	
	
	! This subroutine defines the new aggregate 
	! masses after break-up

	! Input variables
	!-----------------
	! p - pointer for aggregate

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, INTENT(INOUT) :: harvest
   	DOUBLE PRECISION :: x, new, old, daughters, probability
   	DOUBLE PRECISION :: Nn, new_Nn, new_n, extra, add, old_Nn, old_p, new_p
   	DOUBLE PRECISION :: before, after, TEP_t
   	INTEGER :: i

! Initial number of agreagates that an aggregate berak into 
        daughters = 2
        CALL random_number(harvest)
         ! Anusha Commented this
         ! Consturct a power law for the number of fragmants that a particle
         ! breaks into based on Goldthwait et al (2004) as mentioned in Jokulsdottir (2016)
        probability = 0.91*daughters**(-1.56)
        DO WHILE(probability .lt. harvest)
           daughters = daughters + 1
           probability = probability + 0.91*daughters**(-1.56)
        ENDDO

   	! Get the total amount of TEP beofre breakup
   	TEP_t = 0.0
   	i = 1
   	DO WHILE (i .le. 10)
		TEP_t = TEP_t + p%TEP(i,1)*p%n
      	i = i + 1
   	ENDDO
   
   	before = TEP_t
   	IF (p%Nn .gt. daughters) THEN
    	Nn = p%Nn
      	old_Nn = p%Nn
      	extra = MOD(Nn,daughters)
      	Nn = p%Nn - extra

      	new_Nn = Nn/daughters
      	add = extra*p%n/new_Nn
      	new_n = p%n*daughters + add
      	new_p = new_Nn * new_n

      	p%Nn = new_Nn
      	p%n = new_n
      	p%p = new_p
      	p%af = 0
      
      	i = 1
      	DO WHILE (i .le. 10) 
        	p%orgC(i,1) = p%orgC(i,1) * new_Nn/old_Nn
         	p%TEP(i,1) = p%TEP(i,1) * new_Nn/old_Nn
         	i = i + 1
      	ENDDO
	  
      	i = 1
      	DO WHILE (i .le. 4) 
        	p%mineral(i) = p%mineral(i) * new_Nn/old_Nn
         	i = i + 1
      	ENDDO
	  
	  	!This need to be modified accordingly when different size droplets are added
      	i = 1
      	DO WHILE (i .le. 6)
			p%oil(i,1) = p%oil(i,1) * new_Nn/old_Nn
          	i = i + 1
       	ENDDO
	  	  
       	i = 1
       	DO WHILE (i .le. 6)
 			p%sedi(i,1) = p%sedi(i,1) * new_Nn/old_Nn
          	i = i + 1
        ENDDO
      	CALL radius(p)
   	ENDIF 
   	!Total amount of TEP after the break up
   	TEP_t = 0
   	i = 1
   	DO WHILE (i .le. 10)
    	TEP_t = TEP_t + p%TEP(i,1)*p%n
      	i = i + 1
   	ENDDO
   	after = TEP_t
!IF (new_Nn .le. 1 .or. old_Nn .le. 1) print*, '--------hello-----------', p
!IF (p%Nn .le. 3) print*, '***break***', p%r, p%Nn, p%s
!print*, 'breaking', p%r, p%Nn, daughters
   	IF (ABS(before-after) .gt. before) print*, 'Th - break', before, after

END SUBROUTINE break_p
!========================================================================

SUBROUTINE protists(mic, i)
	! What does this do?
	USE the_info
    USE data_sim
   	IMPLICIT NONE
   
   	101 FORMAT (E15.4E2,1X,E15.4E2,1X,E15.4E2,1X,E15.4E2,1X,E15.4E2)

   	INTEGER, INTENT(IN) :: i
   	TYPE(zooplankton), INTENT(INOUT) :: mic
   	INTEGER :: max_i, x
   	DOUBLE PRECISION :: depth

	! depth [cm] at middle of box
   	depth = i*dz - 0.5*dz               

   	mic%Z = mic%P

#ifdef array_size
   	max_i = UBOUND(mic)
   	CALL check_index(i,max_i,x)
   	IF (x .gt. 0) THEN
   		print*, 'OVERWRITE ARRAY in protists'
      	print*, i, max_i
   	ENDIF
#endif

END SUBROUTINE protists
!========================================================================
SUBROUTINE fecal_pellet(p, harvest)
	
   USE the_info
   USE data_sim
   IMPLICIT NONE

   TYPE(agg_part), INTENT(INOUT) :: p
   DOUBLE PRECISION, INTENT(INOUT) :: harvest
   DOUBLE PRECISION :: test1, test2, oldN, before, after
   INTEGER :: i

   oldN = p%Nn
   test1 = p%orgC(1,1)*p%n
   CALL random_number(harvest)
   IF (p%af .eq. 1) THEN
      IF (p%Nn .eq. 0) p%Nn = 100
   ENDIF

   p%n = p%p/p%Nn
   
   i = 1
   DO WHILE (i .le. 10) 
      p%orgC(i,1) = p%orgC(i,1)*p%Nn/oldN
      p%TEP(i,1) = p%TEP(i,1)*p%Nn/oldN
      i = i + 1
   ENDDO
   
   i = 1
   DO WHILE (i .le. 4) 
      p%mineral(i) = p%mineral(i)*p%Nn/oldN
      i = i + 1
   ENDDO
   
   ! This needs to be checked further depending on, if fecal pellets 
   ! are formed with oil mixed aggregates
   i = 1
   DO WHILE (i .le. 6)
	   p%oil(i,1) = p%oil(i,1)*p%Nn/oldN ! Adding this fixed the higer 
	   									! velocites in the particles
	   i = i + 1
   ENDDO
   
   i = 1
   DO WHILE (i .le. 6)
	   p%sedi(i,1) = p%sedi(i,1)*p%Nn/oldN  					
	   i = i + 1
   ENDDO
									    
   test2 = p%orgC(1,1)*p%n
   CALL radius(p)
   CALL density(p)

   IF (ABS(test1-test2) .gt. 1) THEN
      print*, 'fecal pellet', test1, test2, '!'
      print*, p%Nn, p%n, oldN
   ENDIF

   IF (p%Nn .lt. 1 .or. p%n .lt. 0) THEN
      print*, 'fecal pellet' , p%Nn, p%n
   ENDIF

END SUBROUTINE fecal_pellet
!========================================================================
SUBROUTINE sink_2(p)
	
    USE the_info
  	USE data_sim
	
	! This subroutine caculates the 
	! sinking of the aggergate
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: z
   	INTEGER :: i

   	z = p%z + p%w*timestep
   
   	! Getting the flux curve??
   	i = 1
   	DO WHILE (i .le. 9) 
    	IF (p%z .lt. i*1e3 .and. z .gt. i*1e3) THEN 
         	CALL sediment_trap(p,i)    ![i=1..10]
      	ENDIF
      	i = i + 1
   	ENDDO
   	DO WHILE (i .le. 20) 
    	j = i - 9
      	IF (p%z .lt. j*1e4 .and. z .gt. j*1e4) THEN 
        	CALL sediment_trap(p,i)    ![i=11..20]
      	ENDIF
      	i = i + 1
   	ENDDO
   	j = 3
   	DO WHILE (j .le. 8)
    	IF (p%z .lt. j*5e4 .and. z .gt. j*5e4) THEN 
        	CALL sediment_trap(p,i)     ![i=21..26]
      	ENDIF
      	i = i + 1
      	j = j + 1
   	 ENDDO   
   	! Check for the Courant number??
   	! Check if the displacement is greater than the horizontal cell size (Anusha)
   	IF (p%w*timestep .gt. 1000.0) THEN
		write(*,*) 'p%z, w/100.0*24*3600.0, dt, dz,r',p%z, p%w, timestep, p%w*timestep, p%r/1000.
   	ENDIF
   	! Displacement calculation
   	p%z = p%z + p%w*timestep
 
   	IF (p%z .lt. 0) then
	   ! Place the aggregate on the surface
	   p%z = 0.0	!100.0	!changed to zero by Anusha previ
	   ! Depth layer
	   p%idl = 1
	ENDIF

END SUBROUTINE sink_2

SUBROUTINE sink(p, set_mass_flux, buo_mass_flux)
	
	USE the_info
 	USE data_sim
	IMPLICIT NONE
	
	! This subroutine estimates the vertical 
	! sinking/rising of the aggregate and the buoyant or
	! settling mass fluxes across defined 
	! depth levels
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	! set_mass_flux - 
	! buo_mass_flux - 

  	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, DIMENSION(ndz_flux,n_bins,8), INTENT(INOUT) :: set_mass_flux, buo_mass_flux		
   	DOUBLE PRECISION :: z, disp, min_disp, harvest, checking, z_init
   	INTEGER :: i, j, idl_old, idl_new

   	! Displacement of the aggreagate within the timestep
   	disp = p%w*timestep
        idl_old = p%idl
        z_init = p%z

   	! Minimum distance the aggregate need to displace in order to
   	! move to the next depth level
   	if (p%w .lt. 0.0) then 
	   ! For buoyant aggregates
	   min_disp = p%z - dz*(p%idl-1)
   	else
		min_disp = dz*p%idl - p%z
   	endif

   	! Check the displacement and update the depth layer flag of the aggregate
   	if (abs(disp) .ge. min_disp .and. abs(disp) .lt. (min_disp + dz)) then
		!write(*,*) 'Correct aggregae settling'
	    ! Check if the aggregate is at a settling flux estimation depth
		if (mod(p%idl,n_flux).eq. 0 .and. disp .gt. 0.0) then
			!write(*,*) 'Aggregate settling'
	 		! Get the depth level index for flux estimation
	 	   	j = p%idl/n_flux
			!write(*,*) 'Aggregate settling', j, set_mass_flux(1,1,1)
			call mass_flux(p, set_mass_flux, j)
		! Check if the aggregate is at a buoyant flux estimation depth 
		elseif (mod(p%idl,n_flux) .eq. 0 .and. disp .lt. 0.0 .and. p%idl .ne.1) then
!			write(*,*) 'Aggregate rising', p%oil(1,1), p%w
			j = int(p%idl/n_flux)
			call mass_flux(p, buo_mass_flux, j) 
		endif

		! Update the depth layer flag of the aggregate
!		p%idl = p%idl + int(disp/abs(disp))
	! Check if the aggregate skipped an entire depth layer while seeling or rising	
!	elseif (abs(disp) .ge. (min_disp + dz)) then
!		write(*,*) 'Courant condition is violated', disp, '[cm]'
!		write(*,*) 'idl, z, w, dt, dz, r', p%idl, p%z, p%w/100., timestep, p%w*timestep/100., p%r/1000.
		!STOP
	endif

    ! Update the particle depth
    p%z = p%z + disp
    p%idl = int(CEILING(p%z/dz))
    IF (p%w .gt. 0 .and. p%z .gt. 5) THEN
       checking = CEILING(p%z/dz)
      IF (p%idl .ne. checking) print*, p%z, p%idl, checking, int(checking)
    ENDIF
   	! If the aggragte reaches the surface after rising
	! make its depth 0 and put it on the surface depth layer
    if (p%z .lt. 0) then
 	   ! Place the aggregate on the surface
 	   p%z = 0.0	!100.0	!changed to zero by Anusha previ  
	   CALL random_number(harvest) 
 	   ! Depth layer
 	   p%idl = 1	
   endif

   ! IF the aggregate passed a depth bin, then we call sediment_trap
!   idl_new = p%idl
!   IF (idl_new .gt. idl_old) THEN
!      CALL sediment_trap(p, idl_old, idl_new)
!   ENDIF

!#ifdef sedimenttrap
   !Another sediment trap
   IF (z_init .lt. 2000 .and. p%z .gt. 2000) THEN
      ST(1,1) = ST(1,1) + p%orgC(1,1)*p%n  !organic [molC/m2day]
      ST(2,1) = ST(2,1) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,1) = ST(3,1) + p%mineral(3)*p%n  !opal
      ST(4,1) = ST(4,1) + p%TEP(1,1)*p%n  !TEP
      ST(5,1) = ST(5,1) + p%oil(1,1)*p%n  !oil
      ST(6,1) = ST(6,1) + p%sedi(1,1)*p%n  !sedi
   ELSEIF (z_init .lt. 5000 .and. p%z .gt. 5000) THEN
      ST(1,2) = ST(1,2) + p%orgC(1,1)*p%n  !organic [molC/m2day]
      ST(2,2) = ST(2,2) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,2) = ST(3,2) + p%mineral(3)*p%n  !opal
      ST(4,2) = ST(4,2) + p%TEP(1,1)*p%n  !TEP
      ST(5,2) = ST(5,2) + p%oil(1,1)*p%n  !oil
      ST(6,2) = ST(6,2) + p%sedi(1,1)*p%n  !sedi
   ELSEIF (z_init .lt. 10000 .and. p%z .gt. 10000) THEN
      ST(1,3) = ST(1,3) + p%orgC(1,1)*p%n  !organic
      ST(2,3) = ST(2,3) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,3) = ST(3,3) + p%mineral(3)*p%n  !opal
      ST(4,3) = ST(4,3) + p%TEP(1,1)*p%n  !TEP
      ST(5,3) = ST(5,3) + p%oil(1,1)*p%n  !oil
      ST(6,3) = ST(6,3) + p%sedi(1,1)*p%n  !sedi
   ELSEIF (z_init .lt. 20000 .and. p%z .gt. 20000) THEN
      ST(1,4) = ST(1,4) + p%orgC(1,1)*p%n  !organic
      ST(2,4) = ST(2,4) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,4) = ST(3,4) + p%mineral(3)*p%n  !opal
      ST(4,4) = ST(4,4) + p%TEP(1,1)*p%n  !TEP
      ST(5,4) = ST(5,4) + p%oil(1,1)*p%n  !oil
      ST(6,4) = ST(6,4) + p%sedi(1,1)*p%n  !sedi
   ELSEIF (z_init .lt. 50000 .and. p%z .gt. 50000) THEN
      ST(1,5) = ST(1,5) + p%orgC(1,1)*p%n  !organic
      ST(2,5) = ST(2,5) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,5) = ST(3,5) + p%mineral(3)*p%n  !opal
      ST(4,5) = ST(4,5) + p%TEP(1,1)*p%n  !TEP
      ST(5,5) = ST(5,5) + p%oil(1,1)*p%n  !oil
      ST(6,5) = ST(6,5) + p%sedi(1,1)*p%n  !sedi
   ELSEIF (z_init .lt. 100000 .and. p%z .gt. 100000) THEN
      ST(1,6) = ST(1,6) + p%orgC(1,1)*p%n  !organic
      ST(2,6) = ST(2,6) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,6) = ST(3,6) + p%mineral(3)*p%n  !opal
      ST(4,6) = ST(4,6) + p%TEP(1,1)*p%n  !TEP
      ST(5,6) = ST(5,6) + p%oil(1,1)*p%n  !oil
      ST(6,6) = ST(6,6) + p%sedi(1,1)*p%n  !sedi
   ELSEIF (z_init .lt. 150000 .and. p%z .gt. 150000) THEN
      ST(1,7) = ST(1,7) + p%orgC(1,1)*p%n  !organic
      ST(2,7) = ST(2,7) + (p%mineral(1)+p%mineral(2))*p%n  !inorgani
      ST(3,7) = ST(3,7) + p%mineral(3)*p%n  !opal
      ST(4,7) = ST(4,7) + p%TEP(1,1)*p%n  !TEP
      ST(5,7) = ST(5,7) + p%oil(1,1)*p%n  !oil
      ST(6,7) = ST(6,7) + p%sedi(1,1)*p%n  !sedi
   ENDIF
!#endif

END SUBROUTINE sink
!========================================================================
SUBROUTINE respiration(p, bacZoo, lostOrg)
	
   	USE the_info
	USE data_sim
   	IMPLICIT NONE	
	
	! This subroutine calculte the
	! material lost from an aggreagte after
	! respiration
	
	! Input variables
	!-----------------
	! p - pointer for aggregate

   	TYPE(agg_part), INTENT(INOUT) :: p
   	INTEGER, INTENT(IN) :: bacZoo  
   	DOUBLE PRECISION, INTENT(OUT) :: lostOrg
   	DOUBLE PRECISION, DIMENSION(10) :: coef, coefTEP
   	DOUBLE PRECISION, DIMENSION(6) :: coefOil
   	DOUBLE PRECISION :: temp, lostO, organic, temp_effect, TEPC, lostTEP, Oil, lostOil
   	INTEGER :: i, x

	! Get the depth level i
	!CALL find_depth(p%z, i)
	i = p%idl
	
	! Get the temperature at the depth level in deg C
	temp = Temp_z(i) !* stick_param
       !CALL temperature(temp, p%z)
   
	! total amount of organic carbon and TEP
   	organic = 0
   	TEPC = 0
        Oil = 0
	
   	i = 1
        DO WHILE (i .le. 10)
           organic = organic + p%orgC(i,1)
           TEPC = TEPC + p%TEP(i,1)
           i = i + 1
        ENDDO
        ! Total amount of oil 
        i = 1
        DO WHILE (i .le. 6)
           Oil = Oil + p%oil(i,1)
!print*, Oil, p%oil(i,1), 'This is oil', i, p%id
           i = i + 1
        ENDDO
	! coefficients for bacterial (1) and zooplankton (2) respiration
	! of "regular" organic carbon and TEP
   	i = 1
   	temp_effect = exp(DLOG(2.0D00)*(temp-30.0)/10.0)
   	DO WHILE (i .le. 10)
		IF (bacZoo .eq. 1) THEN
!         	coef(i) = 0.1/p%orgC(i,2) * temp_effect * stick_param
!         	coefTEP(i) = 0.1/p%TEP(i,2) * temp_effect * stick_param
         	coef(i) = 0.1/p%orgC(i,2) * temp_effect * 0.3
         	coefTEP(i) = 0.1/p%TEP(i,2) * temp_effect * 0.3
      	 ELSEIF (bacZoo .eq. 2) THEN
         	IF (p%af .eq. 1) THEN      !Micro
            	coef(i) = (1-(2.0**i/2.0**11))/timestep*0.9*temp_effect
            	coefTEP(i) = (1-(2.0**i/2.0**11))/timestep*0.9*temp_effect
         	ENDIF
         	IF (coef(i)*timestep .gt. 1) THEN
			 	print*, 'SOS respiration', temp_effect, coef(i), i
            	coef(i) = 1.0/timestep
         	ENDIF
         	IF (coefTEP(i)*timestep .gt. 1) THEN
			 	print*, 'SOS - TEP', temp_effect, coef(i), i
            	coefTEP(i) = 1.0/timestep
         	ENDIF
      	ENDIF
      i = i + 1
	 ENDDO

   	lostO = 0d00
   	lostOrg = 0d00
        lostOil = 0d00

	IF (organic .gt. 0) THEN
   		i = 1
   		! Track the Oxygen
   		DO WHILE (i .le. 10) 
    		lostO = p%orgC(i,1) * coef(i)*timestep
      		p%orgC(i,1) = p%orgC(i,1) - lostO
      		inventory(1) = inventory(1) - lostO*p%n
      		lost(1) = lost(1) + lostO*p%n
	  		lost_mass_balnce(1) = lost_mass_balnce(1) + lostO*p%n 	  
      		IF (bacZoo .eq. 1) THEN
        		orgC_b(x) = orgC_b(x) + lostO*p%n
			ENDIF
			i = i + 1
         		lostOrg = lostOrg + lostO
   		ENDDO
	ENDIF

	IF (TEPC .gt. 0) THEN
   		i = 1
   	 	DO WHILE (i .le. 10)
      		lostTEP = p%TEP(i,1) * coefTEP(i)*timestep
      		p%TEP(i,1) = p%TEP(i,1) - lostTEP
      		inventory(1) = inventory(1) - lostTEP*p%n
      		lost(5) = lost(5) + lostTEP*p%n
	  		lost_mass_balnce(6) = lost_mass_balnce(6) + lostTEP*p%n
      		i = i + 1
      		lostOrg = lostOrg + lostTEP
   	 	ENDDO
	ENDIF
!print*, Oil, 'OIL'
        IF (Oil .gt. 0) THEN
           i = 1
           DO WHILE (i .le. 6)
!              coefOil(i) = 0.1/p%oil(i,2) * temp_effect * 0.3 * stick_param
              coefOil(i) = 0.1/p%oil(i,2) * temp_effect * 0.3
              lostOil = p%oil(i,1) * coefOil(i)*timestep
              p%Oil(i,1) = p%Oil(i,1) - lostOil
              lost(6) = lost(6) + lostOil*p%n
              lost_mass_balnce(7) = lost_mass_balnce(7) + lostOil*p%n
!print*, lostOil, 'OIL LOST', Oil, i, lost_mass_balnce(7)
              i = i + 1
           ENDDO
        ENDIF
        IF (lostOrg .gt. organic+TEPC) THEN
          print*, 'RESPIRATION - EMERGENCY', p%id
          print*, lostOrg, organic, coef
        ENDIF

        IF (lostOil .gt. Oil) THEN
          print*, 'RESPIRATION - OIL EMERGENCY', p%id
          print*, lostOil, Oil, coefOil(1)
        ENDIF

        CALL radius(p)

END SUBROUTINE respiration
!========================================================================
SUBROUTINE microzoop(p, zoop, harvest)
	
	USE the_info
	USE data_sim
	
	! This subroutine decides if
	! zooplankton eats the agreagete
	! after an encounter
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
   	TYPE(agg_part), INTENT(INOUT) :: p
   	TYPE(zooplankton), INTENT(INOUT) :: zoop   !molC/m^3
   	DOUBLE PRECISION, INTENT(INOUT) :: harvest
   	DOUBLE PRECISION :: zoo, dv, dv1, v, P_enc, orgCfrac, food, fact
   	DOUBLE PRECISION :: oil_frac, sedi_frac, paraD
   	DOUBLE PRECISION :: P_break, apetite, microZ, depth, depth2
   	INTEGER :: i, flag

    !zoopl. dissol.
	paraD = 1
   	microZ = zoop%Z
   	depth = p%z/100.0
   	flag = 0

   	depth2 = depth/15.0
   	fact = 0.1*(depth2**0.5*exp(-depth2/80.0) - 3.9)
   	P_enc = 10**fact * (-1/log(microZ))
   	IF (-1/log(microZ) .gt. 1) THEN
      	print*, 'ERROR: microzoop:', P_enc, microZ, depth
      	P_enc = 0.9
   	ENDIF

   	CALL random_number(harvest)
   	IF (P_enc .gt. harvest) THEN
		! Gets the comparison parameter to be used to decide 
	   	! if the particle is large enough for break up??
      	P_break = DATAN(p%r/1e4) * 0.1
      	CALL random_number(harvest)
      	harvest2 = harvest
      	CALL random_number(harvest)
	  	! Get the fraction of organic carbon from the total weight
      	CALL orgC_fraction(p, orgCfrac)
	  	! Get the fraction of sedi and oil fracirons from the total weight
	  	CALL sedi_oil_fraction(p, oil_frac, sedi_frac) !ALD added this 03/2017	  
      	apetite = orgCfrac 
      	IF (P_break .gt. harvest .and. p%Nn .gt. 1) THEN
         	CALL break_p(p, harvest)
			! Anusha Commented this and added the next line 03/2017
			! If there is oil assumet that zooplankton do not eat the aggregate
			! ELSEIF (apetite .gt. harvest2) THEN		
		 	ELSEIF (apetite .gt. harvest2 .and. oil_frac .le. 0.0 .and. sedi_frac .le. 0.0) THEN
         	   CALL ingestion(p, 1, paraD, harvest)
      	 ENDIF
   	ENDIF
 
END SUBROUTINE microzoop

!========================================================================
SUBROUTINE ingestion(p, flagFP, paraD, harvest)
   
	USE the_info
	USE data_sim
   	IMPLICIT NONE
	
	! This subroutine estimates the 
	! losses from the aggregate when they are 
	! ingested by the marine animals
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
   
   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, INTENT(INOUT) :: paraD, harvest
   	INTEGER, INTENT(IN) :: flagFP !FP: 1=micro 
   	TYPE(agg_part)  :: p_old
   	DOUBLE PRECISION :: organic, lostO, lostC, lostA, orgLost
   	DOUBLE PRECISION :: calcite, aragonite, orgAssim, calcDiss, aragDiss
   	DOUBLE PRECISION :: TEP, lostTEP, S_a, k_size
   	INTEGER :: i, x

   	p_old = p
   	organic = 0
   	i = 1
   	DO WHILE (i .le. 10)
    	organic = organic + p%orgC(i,1)
      	TEP = TEP + p%TEP(i,1)
      	i = i + 1
   	ENDDO
	 
   	IF (organic .lt. 0) print*, 'ERROR: ingestion:', organic, p%orgC(1,1)
   	IF (TEP .lt. 0) print*, 'ERROR: ingestion(TEP)', TEP, p%TEP(1,1)
   	
	calcite = p%mineral(1)
   	aragonite = p%mineral(2)
   
   	!CALL find_depth(p%z, x)
	x = p%idl

   	calcDiss = 0.02*paraD
   	aragDiss = calcDiss/5

	! zooplankton organic carbon and TEP respiration
   	IF (organic .gt. 0) THEN
		IF (flagFP .eq. 1) THEN     !micro
         	p%af = 1
         	CALL respiration(p, 2, lostO)
         	orgC_G1(x) = orgC_G1(x) + lostO*p%n
      	ENDIF 
	  
 	 	!zooplankton calcite dissolution
      	IF (calcite .gt. 0) THEN
         	lostC = lostO/(organic+TEP)*calcite*calcDiss
         	IF (lostC .gt. calcite) THEN
            	lostC = calcite
         	ENDIF

         	p%mineral(1) = p%mineral(1) - lostC
         	lost(2) = lost(2) + lostC*p%n
		 	lost_mass_balnce(2) = lost_mass_balnce(2) + lostC*p%n 
         	inventory(2) = inventory(2) - lostC*p%n
         	calc_G1(x) = calc_G1(x) + lostC*p%n
      	ENDIF
	  
 	 	!zooplankton aragonite dissolution
      	IF (aragonite .gt. 0) THEN
         	k_size = (p%mineral(2)/(p%mineral(2)+1e-10))
         	S_a = 0.1
         	lostA = lostO/(organic+TEP)*aragonite*aragDiss
         !lostA = lostO*aragDiss * (aragonite/(aragonite+1e-10))&
                 !*aragonite*timestep
         !lostA = k_size*k_caco3_a*S_a*p%mineral(1)*timestep
         IF (lostA .ge. aragonite) THEN
            lostA = aragonite
         ENDIF

         p%mineral(2) = p%mineral(2) - lostA
         lost(2) = lost(2) + lostA*p%n
		 lost_mass_balnce(3) = lost_mass_balnce(3) + lostA*p%n
         inventory(2) = inventory(2) - lostA*p%n
         calc_G1(x) = calc_G1(x) + lostA*p%n
      ENDIF
  
  	ENDIF
   
	CALL fecal_pellet(p, harvest)

	i = 1
	DO WHILE (i .lt. 10) 
   		IF (p%orgC(i,1) .lt. 0) THEN
      	  p%orgC(i,1) = 0.0
   		ENDIF
   	IF (p%orgC(i,1) .lt. 1d-50 .and. p%orgC(i,1) .gt. 0.0) THEN
    	inventory(1) = inventory(1) - p%orgC(i,1)*p%n 
      	lost(1) = lost(1) + p%orgC(i,1)*p%n
	  	lost_mass_balnce(1) = lost_mass_balnce(1) + p%orgC(i,1)*p%n	  
      	p%orgC(i,1) = 0.0
   	ENDIF
   	i = i + 1
	ENDDO

END SUBROUTINE ingestion

!========================================================================
!========================================================================
SUBROUTINE photolysis(p)
   	
	USE the_info
	USE data_sim
	
	! Need to figure out what is 
	! happening here
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	!Dissolution of Organic and minerals

   	TYPE(agg_part), INTENT(INOUT) :: p
	DOUBLE PRECISION :: rho_sw_m
		
	!Get the density of sea water at the aggregate depth level 
   	rho_sw_m = Dens_z(p%idl)
	
   	IF (p%rho .lt. rho_sw_m) THEN
    	p%r = 0.1
      	CALL dissolve(p)
   	ENDIF
   
END SUBROUTINE photolysis
!========================================================================
 SUBROUTINE temperature(temp, depth)
    USE the_info
    USE  data_sim
 ! depth is in cm, temp in C

    DOUBLE PRECISION, INTENT(OUT) :: temp
    DOUBLE PRECISION, INTENT(IN) :: depth
!    DOUBLE PRECISION, DIMENSION(2,37) :: iTemp   !eqPac
    DOUBLE PRECISION, DIMENSION(2,3) :: iTemp   !NA or SO
    DOUBLE PRECISION :: z
    INTEGER :: i

    iTemp(1,:) = (/ paraT,4.0D00,2.0D00 /)
    iTemp(2,:) = (/ 0,1000,4000 /)

    z = depth/100.0
    i = 1
    DO WHILE ( z .ge. iTEMP(2,i) )
       i = i + 1
    ENDDO

    CALL interpolate(iTemp(1,i-1),iTemp(1,i),iTemp(2,i-1),iTemp(2,i),&
         temp,z)

    IF (temp .gt. 35) THEN
       print*, 'temperature greater than 35'
       print*, iTEMP(1,1), iTemp(2,1), temp
    ENDIF

 END SUBROUTINE temperature
!========================================================================
SUBROUTINE Temperature_Measured(temp, z)
	
	USE the_info
	USE data_sim
	 
	! Get the temperature at the defined depth
	! interpolated from measured data profile
	! Input Parameters
	!------------------
	! T : float
	!	Temperature (deg C)
	! z : float
	!	Depth (cm)
	
   	DOUBLE PRECISION, INTENT(OUT) :: temp
   	DOUBLE PRECISION, INTENT(IN) :: z
   	INTEGER :: i

	! Get the depth level i
	CALL find_depth(z, i)
	
	! Get the temperature
	temp = Temp_z(i) !* stick_param


END SUBROUTINE Temperature_Measured
!========================================================================
SUBROUTINE aging(p)
	
   	USE the_info
    USE data_sim
	
	! Need to figure out how to define the amount of organic carbon
	! is old
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Aging of oil needs to be added to this 

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION :: older
   	INTEGER :: i

   	i = 9
   	older = 0

   	DO WHILE (i .ge. 1) 
    	IF (i .eq. 1) THEN
        	older = p%orgC(i,1) * timestep/(p%orgC(i,2))
      	ELSE
        	older = p%orgC(i,1) * timestep/(p%orgC(i,2)-p%orgC(i-1,2))
      	ENDIF
      	IF (older .gt. p%orgC(i,1)) THEN
        	p%orgC(i+1,1) = p%orgC(i+1,1) + p%orgC(i,1)
         	p%orgC(i,1) = 0.0
      	ELSE
        	p%orgC(i,1) = p%orgC(i,1) - older
         	p%orgC(i+1,1) = p%orgC(i+1,1) + older
      	ENDIF
      	i = i - 1
   	 ENDDO
   
   	i = 9
   	older = 0
   	DO WHILE (i .ge. 1) 
    	IF (i .eq. 1) THEN
        	older = p%TEP(i,1) * timestep/(p%TEP(i,2))
      	ELSE
        	older = p%TEP(i,1) * timestep/(p%TEP(i,2)-p%TEP(i-1,2))
      	ENDIF
      	IF (older .gt. p%TEP(i,1)) THEN
        	p%TEP(i+1,1) = p%TEP(i+1,1) + p%TEP(i,1)
         	p%TEP(i,1) = 0.0
      	ELSE
        	p%TEP(i,1) = p%TEP(i,1) - older
         	p%TEP(i+1,1) = p%TEP(i+1,1) + older
      	ENDIF
      	i = i - 1
   	ENDDO

END SUBROUTINE aging

!========================================================================
SUBROUTINE orgC_fraction(p, dwOrgCfrac)

	USE the_info
  	USE data_sim	

	! This subroutine calculates the organic carbon dry weight as a fraction
	! of total weight
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Output variables
	!-----------------
	! dwOrgCfrac - dry weigt of organic carbon

   	TYPE(agg_part), INTENT(INOUT) :: p
   	DOUBLE PRECISION, INTENT(INOUT) :: dwOrgCfrac
   	DOUBLE PRECISION :: organic, dw_mineral, dw_orgC, dw_oil, dw_sedi, dw_total
   	DOUBLE PRECISION ::  TEP, dw_TEP

	! amount of organic carbon
   	organic = 0.0
	! amount of TEP
   	TEP = 0.0
	! weight of oil
   	dw_oil = 0.0
	! dry weight of sediments
   	dw_sedi = 0.0
	! weight of mineral
   	dw_mineral = (p%mineral(1)+p%mineral(2))*mw_co + &
             p%mineral(3)*mw_si + p%mineral(4)*mw_li
			 
   	i = 1
   	DO WHILE (i .le. 10) 
    	organic = organic + p%orgC(i,1)
      	TEP = TEP + p%TEP(i,1)
      	i = i + 1
   	ENDDO
   
   	!dry weight of organic carbon
   	dw_orgC = organic*mw_orgC/OMCfrac
	! weight of TEP
   	dw_TEP = TEP*mw_TEP
   
   	i = 1
   	DO WHILE (i .le. 6)
		dw_oil = dw_oil + p%oil(i,1)
       	i = i + 1
    ENDDO	   
  
    i = 1
    DO WHILE (i .le. 6)
 		dw_sedi = dw_sedi + p%sedi(i,1)
        i = i + 1
    ENDDO	   

	! Dry weight of total aggregate
   	dw_total = dw_orgC + dw_mineral + dw_TEP + dw_oil + dw_sedi
	! Get the fraction of dry weight of irganic carbon in the aggregate
   	dwOrgCfrac = dw_orgC / dw_total

END SUBROUTINE orgC_fraction

!========================================================================
SUBROUTINE sedi_oil_fraction(p, oil_frac, sedi_frac)
    	
	USE the_info
  	USE data_sim	
	
	! This subroutine calculates the oil dry weight and sedi dry weight 
	! as a fraction of total weight

	! Input variables
	!-----------------
	! p - pointer for aggregate

	! Output variables
	!-----------------
	! p%oil_frac - mass fraction of oil
	! oil_frac - mass fraction of sediment
	
   	TYPE(agg_part), INTENT(INOUT) :: p
  	DOUBLE PRECISION, INTENT(INOUT) :: oil_frac, sedi_frac
   	DOUBLE PRECISION :: organic, dw_mineral, dw_orgC, dw_oil, dw_sedi, dw_total
   	DOUBLE PRECISION :: TEP, dw_TEP

	! amount of oraganic carbon
   	organic = 0.0
   	! amount of TEP
	TEP = 0.0
	!weight of oil
   	dw_oil = 0.0
	! dry weight of sediment
   	dw_sedi = 0.0
	! weight of minerals
   	dw_mineral = (p%mineral(1)+p%mineral(2))*mw_co + &
             p%mineral(3)*mw_si + p%mineral(4)*mw_li
			 
	i = 1
   	DO WHILE (i .le. 10) 
    	organic = organic + p%orgC(i,1)
      	TEP = TEP + p%TEP(i,1)
      	i = i + 1
   	 ENDDO
   
   	dw_orgC = organic*mw_orgC/OMCfrac
   	dw_TEP = TEP*mw_TEP
   
   	i = 1
   	DO WHILE (i .le. 6)
		dw_oil = dw_oil + p%oil(i,1)
       	i = i + 1
    ENDDO	   
  
    i = 1
    DO WHILE (i .le. 6)
 		dw_sedi = dw_sedi + p%sedi(i,1)
        i = i + 1
     ENDDO	   

	! Total dry weight of the aggregate
   	dw_total = dw_orgC + dw_mineral + dw_TEP + dw_oil + dw_sedi

   	! Get the mass fraction of oil
	oil_frac = dw_oil / dw_total
	! Get the mass fraction of sediment
   	sedi_frac = dw_sedi / dw_total

END SUBROUTINE sedi_oil_fraction

!========================================================================
SUBROUTINE dissolution(p)

 	USE the_info
	USE data_sim
   	IMPLICIT NONE	
	
	! This subroutine calculates the disslution of
	! different marine snow components
	
	! Input variables
	!-----------------
	! p - pointer for aggregate
	
	! Oil dissolution needs to be added yet


   	TYPE(agg_part), INTENT(INOUT) :: p
   	INTEGER :: x, i
   	DOUBLE PRECISION :: R, Rd, Ta, temp, S_c, S_a
   	DOUBLE PRECISION :: lostC, lostS, lostA, k_size
   	DOUBLE PRECISION :: mass_oil, voil, mass_sedi, vsedi
   	DOUBLE PRECISION :: calc, clay, opal, OrgC, tvol
   	DOUBLE PRECISION :: vcalc, vclay, vopal, vorgC
   	DOUBLE PRECISION :: k_depth, max_z
   
   	lostA = 0.0
   	lostC = 0.0
   	lostS = 0.0
   	S_c = 0.0
   	S_a = 0.0

   	IF (p%mineral(1)+p%mineral(2)+p%mineral(3) .gt. 0) THEN
		! Get the depth level
		!CALL find_depth(p%z, x)
		x = p%idl
		
		! Get the tmperature withn the depth level i deg C
		temp = Temp_z(x) !* stick_param
		
      	CALL deltaCO3(p%z, S_c, S_a)

		! Convert the tepmpearture to Kelvin
      	Ta = temp+273.15
      	R =  1.32e16*exp(-11481/Ta)  ![day^-1]
      	Rd = R/86400                 ![s^-1]

		! Commented by Anusha 04/2017
      	!CALL find_depth(p%z, x)

		! CaCO3 thermodynamics
 	   	! if amount of CaCO3 in agg less than in 1/2 coccolithoph
 	  	! then the whole thing dissolves, otherwise according to thermodyn.
      	IF (S_c .gt. 0) THEN
         	IF (p%mineral(1) .gt. 2.5e-12) THEN
            	k_size = (p%mineral(1)/(p%mineral(1)+1e-10))
            	lostC = k_caco3_c*S_c*p%mineral(1)*timestep

            	IF (p%mineral(1) .lt. 2.5e-12 .or. p%mineral(1) .lt. lostC) THEN  
               	 lostC = p%mineral(1)
            	ENDIF

            	p%mineral(1) = p%mineral(1) - lostC
            	inventory(2) = inventory(2) - lostC*p%n
            	lost(2) = lost(2) + lostC*p%n
				lost_mass_balnce(2) = lost_mass_balnce(2) + lostC*p%n 			
            	calcC_t(x) = calcC_t(x) + lostC*p%n
         	ENDIF
	 	ENDIF
	  
	 	! aragonite dissolution       
     	IF (S_a .gt. 0) THEN
		 	IF (p%mineral(2) .gt. 2.5e-12) THEN
         		k_size = (p%mineral(2)/(p%mineral(2)+1e-10))
            	lostA = k_caco3_a*S_a*p%mineral(2)*timestep
            	IF (lostA .gt. p%mineral(2)) THEN
            		print*, p%mineral(2), lostA, S_a
               		print*, '-------mineral(2)-----------', p%z/100
            	ENDIF
           	 	IF (p%mineral(2) .lt. 2.5e-12 .or. lostA .gt. p%mineral(2)) THEN  
               		lostA = p%mineral(2)
            	ENDIF

            	p%mineral(2) = p%mineral(2) - lostA
            	inventory(2) = inventory(2) - lostA*p%n
            	lost(2) = lost(2) + lostA*p%n
		lost_mass_balnce(3) = lost_mass_balnce(3) + lostA*p%n 
            	calcA_t(x) = calcA_t(x) + lostA*p%n
        	 ENDIF
      	ENDIF

          ! zooplankton CaCO3 dissolution
          k_depth = 0
          max_z = 100000 !cm
          IF (p%mineral(1) .gt. 0) THEN
             IF (p%z .le. max_z) THEN
!                k_depth = stick_param * sin(3.1415*p%z/max_z) / timestep
                k_depth = 0.002 * sin(3.1415*p%z/max_z) / timestep
             ELSEIF (p%z .gt. max_z) THEN
                k_depth = 0
             ENDIF
             lostC = k_depth*p%mineral(1) * timestep
             p%mineral(1) = p%mineral(1) - lostC
             inventory(2) = inventory(2) - lostC*p%n
             lost(2) = lost(2) + lostC*p%n
             lost_mass_balnce(2) = lost_mass_balnce(2) + lostC*p%n 
             IF (p%mineral(1) .lt. 2.5e-12 .or. p%mineral(1) .lt. lostC) THEN  
                lostC = p%mineral(1)
             ENDIF
          ENDIF

          ! zooplankton CaCO3_aragonite dissolution
          k_depth = 0
          max_z = 100000 !cm
          IF (p%mineral(2) .gt. 0) THEN
             IF (p%z .le. max_z) THEN
!                k_depth = stick_param * sin(3.1415*p%z/max_z) / timestep
                k_depth = 0.002 * sin(3.1415*p%z/max_z) / timestep
             ELSEIF (p%z .gt. max_z) THEN
                k_depth = 0
             ENDIF
             lostA = k_depth*p%mineral(2) * timestep
             p%mineral(2) = p%mineral(2) - lostA
             inventory(2) = inventory(2) - lostA*p%n
             lost(2) = lost(2) + lostA*p%n
             lost_mass_balnce(3) = lost_mass_balnce(3) + lostA*p%n 
             IF (p%mineral(1) .lt. 2.5e-12 .or. p%mineral(1) .lt. lostC) THEN  
                lostC = p%mineral(1)
             ENDIF
          ENDIF




	  	! opal dissolution
      	IF (p%mineral(3) .gt. 0) THEN
!      		k_size = (p%mineral(3)/(p%mineral(3)+stick_param))
      		k_size = (p%mineral(3)/(p%mineral(3)+2e-10))
        	lostS = k_size*p%mineral(3)*Rd*timestep
        	p%mineral(3) = p%mineral(3) - lostS
        	inventory(3) = inventory(3) - lostS*p%n
        	lost(3) = lost(3) + lostS*p%n
		lost_mass_balnce(4) = lost_mass_balnce(4) + lostS*p%n 		 
        	opal_t(x) = opal_t(x) + lostS*p%n
      	ENDIF

	  	! volume of solid stuff
      	i = 1
      	orgC = 0d00
      	DO WHILE (i .le. 10) 
      		orgC = orgC + p%orgC(i,1) 
        	i = i + 1
      	  ENDDO
      	calc = p%mineral(1) + p%mineral(2)
      	opal = p%mineral(3)
      	clay = p%mineral(4)
	  
      	i = 1
      	mass_oil = 0.0	
      	DO WHILE (i .le. 6) 	    
	  		mass_oil = mass_oil + p%oil(i,1) 
        	i = i + 1
     	ENDDO
  
    	i = 1
    	mass_sedi = 0.0	
    	DO WHILE (i .le. 6) 	    
  			mass_sedi = mass_sedi + p%sedi(i,1) 
       		i = i + 1
    	ENDDO
	   
		! change moles to volume
		vorgC = (orgC * mw_orgC / OMCfrac) / rho_orgC
    	vcalc = calc * mw_co / rho_co
    	vopal = opal * mw_si / rho_si
    	vclay = clay * mw_li / rho_li
		voil = mass_oil / rho_oil
		vsedi = mass_sedi / rho_sedi

		! Get the total volume of the aggregate
		tvol = vorgC + vcalc + vopal + vclay + voil + vsedi
   
		IF (lostS+lostC+lostA .gt. 0) THEN
			!write(*,*) 'dissolution'
			!write(*,*) 'p%r bofore',p%r
        	CALL radius(p)
			!write(*,*) 'p%r after',p%r
		ENDIF
   	
	ENDIF

	IF (p%mineral(1) .lt. 0) print*, 'dissol', p%mineral(2), p%TEP(1,1)

END SUBROUTINE dissolution
!=======================================================================

SUBROUTINE deltaCO3(depth, S_c, S_a)
	
	USE the_info
	USE data_sim

   	DOUBLE PRECISION, INTENT(IN) :: depth
   	DOUBLE PRECISION, INTENT(INOUT) :: S_c, S_a
   	DOUBLE PRECISION :: co3ion, co3ion1, co3ion2, co3ion3, co3iond
   	DOUBLE PRECISION :: co3satC, co3satA, z, paraC

   	z = depth/100  !m
	
	! CaCO3
	paraC = 220 
	
   	co3ion1 = paraC
   	co3ion2 = paraC - 140
   	co3ion3 = paraC - 140 !paraC/5 + 40
   	co3iond = 1000
!   co3ion1 = 210 !PACIFIC (umol/kg)
!   co3ion2 = 65 
!   co3ion3 = 80
!   co3iond = 1000

!   co3ion1 = 250 !NATLANTIC(umol/kg)
!   co3ion2 = 150 
!   co3ion3 = 100
!   co3iond = 1000

!   co3ion1 = 230 !SOUTHERN(umol/kg)
!   co3ion2 = 100 
!   co3ion3 = 90
!   co3iond = 1000

   	IF (z .lt. 2500) THEN
      	co3satC = 0.01104*z + 41.8
      	co3satA = 0.0168*z + 65
   	ELSE
      	co3satC = 0.01692*z + 27.1
      	co3satA = 0.024*(z - 2500) + 107
   	ENDIF

   	IF (z .lt. co3iond) THEN
      	slope = (z - co3iond)/((-co3iond)*2/pi)
      	co3ion = co3ion1 + (co3ion2-co3ion1)*cos(slope)
   	ELSE
      	co3ion = co3ion2+(co3ion3-co3ion2)*tan(0.86*(z-co3iond)/(seafloor/100))
   	ENDIF
   
   	IF (co3ion .gt. co3satC) THEN
      	S_c = 0
   	ELSE
      	S_c = (1 - co3ion/co3satC)**eta_c
   	ENDIF

   	IF (co3ion .gt. co3satA) THEN
      	S_a = 0
   	ELSE
      	S_a = (1 - co3ion/co3satA)**eta_a
   	ENDIF

END SUBROUTINE deltaCO3

!========================================================================
SUBROUTINE viscosity(depth, nu)

	USE the_info
	USE data_sim	
	
	! Dynamic viscosity in units g/cms
	! Anusha: I think this is in kg/ms ?? Check!!!

	DOUBLE PRECISION, INTENT(IN) :: depth
   	DOUBLE PRECISION, INTENT(OUT) :: nu
   	DOUBLE PRECISION :: temp
	INTEGER :: i
   
	! Get the depth level i
	CALL find_depth(depth, i)

	! Get the temperature at the depth level in deg C
	temp = Temp_z(i) !* stick_param
	
   	! Get the temperature within the depth level
   	CALL Temperature_Measured(temp, depth)
!        CALL temperature(temp, depth)

   	nu = 5e-6 * exp(2250/(temp+273.15))
   	IF (nu .lt. 1e-3 .or. nu .gt. 1) THEN
    	print*, 'viscosity'
      	print*, depth, temp, nu
   	ENDIF

END SUBROUTINE


SUBROUTINE seawater_density(density)
	
	! Have function in tamoc that can be used directly
	DOUBLE PRECISION :: density
	
	density = 1026
	
END SUBROUTINE

SUBROUTINE Energy_dissipation(depth, epsil)
	
	! Calculates the energy dissipation epsil
	! at a defined depth and a surface wind 
	! based on MacKenzie and Leggett (1993)
	! Added by Anusha 05/2017
	
	! Input variables
	!----------------
	! depth - location depth (m)
	
	!Output variables
	!----------------
	!epsil - energy dissipation (W/kg or m2/s3)
	
	USE data_stations
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN) :: depth	
	DOUBLE PRECISION, INTENT(OUT) :: epsil 

	! Calculate the energy dissipation in m2/s3
	epsil = 5.82e-9*U_wind**3*(depth)**(-1)

END SUBROUTINE	

SUBROUTINE Average_shear_rate(kin_vis, depth, ave_gamma)

	! Average Shear rate is  estimated based on 
	! Jakson (2001) and Pruppacher and Klett (1980)
	! Added by Anusha 05/2017	
	! I used the depht of each layer in the formulaiton instead 
	! of the mixing layer dpeth suggested in the above papers.
	! Seems its is ok for 300 m depth simulation
	
	! Input variables
	!-----------------
	! kin_vis - Kinemetic viscosity (m^2/s)
	
	! Output variables
	!-----------------
	! Ave_gamma - Average shear rate (s-1)
	
	USE data_stations
	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN) :: kin_vis, depth
	DOUBLE PRECISION, INTENT(OUT) :: ave_gamma
	
	! Calculate the average shear rate
	ave_gamma = 7.62e-5*U_wind**1.5*kin_vis**(-0.5)*depth**(-0.5)
	!ave_gamma = 1.53e-4*U_wind**1.5*kin_vis**(-0.5)*Z_mld**(-0.5)	
	!write(*,*) 'ave_gamma',ave_gamma, kin_vis, U_wind, Z_mld

END SUBROUTINE


SUBROUTINE Kolmogorov_length_scale(kin_vis, epsil, eta)
	
	USE data_stations
	IMPLICIT NONE
	
	! Calculates the Kolmogorov length scale when 
	! the Kinemetic viscosity and the energy dissipation rate
	! is given
	
	! Input variables
	!-----------------
	! kin_vis - Kinemetic viscosity (m^2/s)
	! epsil - energy dissipation (W/kg or m2/s3)
	
	! Output variables
	!-----------------
	! eta - Kolmogorov length scale (m)	
	
	DOUBLE PRECISION, INTENT(IN) :: kin_vis, epsil
	DOUBLE PRECISION, INTENT(OUT) :: eta
	
	! Calcualte the Kolmogorov length scale in m
	eta = (kin_vis**3/epsil)**(0.25)
	
END SUBROUTINE

SUBROUTINE update_time(p) 
	
	USE the_info
	USE  data_sim	
	IMPLICIT NONE
   	TYPE(agg_part), INTENT(INOUT) :: p	
	
	! Udate the time of the aggregate age in the water 
	
	!Input 
	!-----
	
	!Output
	!----
	
	!Updatae the time stamp of the aggregate
	p%t = p%t + timestep*1.0
	
END SUBROUTINE



	
