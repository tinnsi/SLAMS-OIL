! """
! Equations of State for Seawater
! ===============================
!
! This module provides a few simple equations of state that can be used for
! seawater.  Some are very good (e.g., density) and others are just taken from
! tables in the Handbook of Physics and Chemistry for pure water (e.g., surface
! tension).  Eventually, these should be replaced with routines from the
! official seawater equation of state.
!
! """
! Functions are extracted from Texas A&M Oil Spill Calculator (TAMOC) and - Converted to fortran functions.
! A.Dissanayake, March 2017, University of Georgia 

SUBROUTINE density_sw(T_in, S_in, P_in, rho_sw)
!     """
!     Computes the density of seawater from Gill (1982)
!
!     Computes the density of seawater using the equation of state in Gill
!     (1982), *Ocean-Atmosphere Dynamics*, Academic Press, New York.  The
!     equations for this code are taken from Appendix B in Crounse (2000).
!
!     Parameters
!     ----------
!     T : float
!         temperature (K)
!     S : float
!         salinity (psu)
!     P : float
!         pressure (Pa)
!
!     Returns
!     -------
!     rho : float
!         seawater density (kg/m^3)
!
!     """

   IMPLICIT NONE
   DOUBLE PRECISION, INTENT(IN) :: T_in, S_in, P_in    
   DOUBLE PRECISION :: T, S, P   
   DOUBLE PRECISION, INTENT(OUT) :: rho_sw
   DOUBLE PRECISION :: g, K
   DOUBLE PRECISION :: rho_sw_0
   
	! Acceleration of gravity (m/s^2)
   g = 9.81     
     
   ! Convert T to dec C and P to bar
    T = T_in - 273.15
    P = P_in * 1.e-5
	
	! Salinity
	S = S_in
    
    ! Compute the density at atmospheric pressure
    rho_sw_0 = (999.842594 + 6.793952e-2 * T - 9.095290e-3 * T**2 &
                + 1.001685e-4 * T**3 - 1.120083e-6 * T**4 + 6.536332e-9 * T**5 &
                + 8.24493e-1 * S - 5.72466e-3 * S**(3./2.) + 4.8314e-4 * S**2 &
                - 4.0899e-3 * T*S + 7.6438e-5 * T**2 * S - 8.2467e-7 * T**3 * &
                S + 5.3875e-9 * T**4 * S + 1.0227e-4 * T * S**(3./2.) &
                - 1.6546e-6 * T**2 * S**(3./2.))
    
     ! Compute the pressure correction coefficient
    K = (19652.21 + 148.4206 * T - 2.327105 * T**2 + 1.360477e-2 * T**3 &
         - 5.155288e-5 * T**4 + 3.239908 * P + 1.43713e-3 * T * P &
         + 1.16092e-4 * T**2 * P - 5.77905e-7 * T**3 * P &
         + 8.50935e-5 * P**2 - 6.12293e-6 * T * P**2 &
         + 5.2787e-8 * T**2 * P**2 + 54.6746 * S - 0.603459 * T * S &
         + 1.09987e-2 * T**2 * S - 6.1670e-5 * T**3 * S &
         + 7.944e-2 * S**(3./2.) + 1.64833e-2 * T * S**(3./2.) &
         - 5.3009e-4 * T**2 * S**(3./2.) + 2.2838e-3 * P * S &
         - 1.0981e-5 * T * P * S - 1.6078e-6 * T**2 * P * S &
         + 1.91075e-4 * P * S**(3./2.) - 9.9348e-7 * P**2 * S &
         + 2.0816e-8 * T * P**2 * S + 9.1697e-10 * T**2 * P**2 * S)
		 
	! Compute the seawater density		 
	rho_sw = rho_sw_0 / (1.0 - P / K)
    
END SUBROUTINE 

SUBROUTINE  dynamic_viscosity_sw(Ta, Sa, Pa, mu)
!     """
!     Compute the viscosity of seawater
!
!     Evaluates the viscosity of seawater as a function of temperature,
!     salinity, and pressure.  The equation accounting for temperature and
!     salinity is from Sharqawy et al. (2010).  The pressure correction is
!     from McCain (1990), equation B-75 on page 528.
!
!     Parameters
!     ----------
!     T : float
!         temperature (K)
!     S : float
!         salinity (psu)
!     P : float
!         pressure (Pa)
!
!     Returns
!     -------
!     mu : float
!         dynamic viscosity of seawater (Pa s)
!
!     """

   	IMPLICIT NONE
	
	DOUBLE PRECISION, INTENT(IN) :: Ta, Sa, Pa
	DOUBLE PRECISION :: Sc, Pc, Tc
   	DOUBLE PRECISION, INTENT(OUT) :: mu	
   	DOUBLE PRECISION :: g, A, B
   	DOUBLE PRECISION :: rho_sw, rho_sw_0, mu_w, mu_0
	DOUBLE PRECISION, DIMENSION(10) :: a_coef
	
	! Acceleration of gravity (m/s^2)
   	g = 9.81  
   
    ! The following equations use Temperature in deg C
    Tc = Ta - 273.15
    
    ! Get the fit coefficients
    a_coef(1:10) = (/1.5700386464E-01, 6.4992620050E+01, -9.1296496657E+01, &
                  4.2844324477E-05, 1.5409136040E+00, 1.9981117208E-02, &
                  -9.5203865864E-05, 7.9739318223E+00, -7.5614568881E-02, &
                  4.7237011074E-04/)
                  
    ! Compute the viscosity of pure water at given temperature
    mu_w = a_coef(4) + 1./(a_coef(1) * (Tc + a_coef(2))**2 + a_coef(3))
    
    ! Correct for salinity
    Sc = Sa / 1000.
    A = a_coef(5) + a_coef(6) * Tc + a_coef(7) * Tc**2
    B = a_coef(8) + a_coef(9) * Tc + a_coef(10) * Tc**2
    mu_0 = mu_w * (1. + A * Sc + B * Sc**2)
    
    ! And finally for pressure
    Pc = Pa * 0.00014503773800721815
    mu = mu_0 * (0.9994 + 4.0295e-5 * Pc + 3.1062e-9 * Pc**2)
    
END SUBROUTINE

SUBROUTINE  dynamic_viscosity_profile(Ta, Sa, Pa, dyna_vis, n_dz)
!     """
!     Compute the viscosity of seawater
!
!     Evaluates the viscosity of seawater as a function of temperature,
!     salinity, and pressure.  The equation accounting for temperature and
!     salinity is from Sharqawy et al. (2010).  The pressure correction is
!     from McCain (1990), equation B-75 on page 528.
!
!     Parameters
!     ----------
!     T : float
!         temperature (K)
!     S : float
!         salinity (psu)
!     P : float
!         pressure (Pa)
!
!     Returns
!     -------
!     mu : float
!         dynamic viscosity of seawater (Pa s)
!
!     """

	USE data_stations
	
   	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(n_dz), INTENT(INOUT) :: dyna_vis
	DOUBLE PRECISION, DIMENSION(n_dz), INTENT(IN) :: Sa, Ta, Pa
   	DOUBLE PRECISION :: mu, g
	INTEGER :: i, n_dz

!     ! Define an array for storing the pressure
! 	allocate(dyna_vis(n_dz), STAT = AllocateStatus)
! 	if (AllocateStatus /= 0) stop "Not enough memory"

	! Acceleration of gravity (m/s^2)
   	g = 9.81
! 	write(*,*) 'g', g

   	do i = 1, n_dz
		CALL dynamic_viscosity_sw(Ta(i), Sa(i), Pa(i), mu)
		dyna_vis(i) = mu
	end do

    write(*,*) 'End of the Dynamic viscosity calcualtion'

END SUBROUTINE

SUBROUTINE interfacial_tension(T, S, sigma)
!     """
!     Compute the surface tension of seawater
!
!     Evaluates the surface tension of seawater as a function of temperature and
!     salinity following equations in Sharqawy et al. (2010), Table 6.
!
!     Parameters
!     ----------
!     T : float
!         temperature (K)
!     S : float
!         salinity (psu)
!
!     Returns
!     -------
!     sigma : float
!         interfacial tension of air in seawater (N/m)
!
!     """

   IMPLICIT NONE
   DOUBLE PRECISION :: T, S
   DOUBLE PRECISION, INTENT(OUT) :: sigma
   DOUBLE PRECISION :: sigma_w
   
    ! Equations in Sharqawy using deg C and g/kg
    T = T - 273.15
    S = S / 1000.
    
    ! Use equations (27) for pure water surface tension (N/m)
    sigma_w = 0.2358 * (1. - (T + 273.15) / 647.096)**1.256 * (1. - 0.625 * &
              (1. - (T + 273.15) / 647.096))
    
    ! Equation (28) gives the salinity correction
    if (T < 40) then
        ! Salinity correction only valid [0, 40] deg C
        sigma = sigma_w * (1. + (0.000226 * T + 0.00946) * log(1. + 0.0331 * &
                S))
    else
        ! No available salinity correction for hot cases
        sigma = sigma_w
		
	end if
    
END SUBROUTINE

SUBROUTINE conductivity_sw(T, S, P, k_sw)
!     """
!     Compute the thermal conductivity of seawater
!
!     Evaluates the thermal conductivity of seawateras a function of
!     temperature, pressure, and salinity following equations in Sharqawy et
!     al. (2010), Table 4.
!
!     Parameters
!     ----------
!     T : float
!         temperature (K)
!     S : float
!         salinity (psu)
!     P : float
!         pressure (Pa)
!
!     Returns
!     -------
!     k : float
!         thermal conductivity of seawater (W/(mK))
!
!     Notes
!     -----
!     Table 4 provides three equations.  Equation (14) is valid for temperatures
!     up to 60 deg C, but I could not reproduce the values in the paper.  Hence,
!     the slightly less-accurate equation (13) is used for temperatures above
!     30 deg C.
!
!     """
	
	IMPLICIT NONE
	
	DOUBLE PRECISION :: T, S, P
	DOUBLE PRECISION :: T_68
	DOUBLE PRECISION, INTENT(OUT) :: k_sw
	
    ! Thermal conductivity equations use T_68 in deg C
    T_68 = (T - 0.0682875) / (1.0 - 0.00025)
    T_68  = T_68 - 273.15
    
    ! Salinity is in g/kg and pressure in MPa
    S = S / 1000.
    P = P * 1e-6
    
    ! Compute the thermal conductivity from Table 4
    if (T_68 < 30.) then
        ! Equation (15)
        k_sw = 0.55286 + 3.4025e-4 * P + 1.8364e-3 * T_68 - 3.3058e-7 * &
               T_68**3
    else
        ! Equation (13)
        k_sw = 10.**(log10(240. + 0.0002 * S) + 0.434 * (2.3 - (343.5 + &
               0.037 * S) / (T_68 + 273.15)) * (1. - (T_68 + 273.15)/ &
               (647. + 0.03 * S)) ** 0.333) / 1000.
	endif
	
END SUBROUTINE
    
SUBROUTINE heat_capacity_sw(cp)
!     """
!     Compute the heat capacity of seawater at fixed conditions
!
!     Per Figure 5 in Sharqawy et al. (2010), the heat capacity of seawater
!     only varies +/- 5 percent over practical temperatures and salinities
!     for deepwater blowouts.  If we let heat capacity depend on temperature
!     and or salinity, computing the temperature of water given the total
!     heat becomes an implicit calculation.  This is a problem for the plume
!     models.  As a result, we choose to set the heat capacity to that of
!     seawater at 10 deg C and 34.5 psu.
!
!     Returns
!     -------
!     cp : float
!         heat capacity of seawater (J/(kg K))
!
!     Note
!     ----
!     This approximation is valid since we have treated cp to be a constant in
!     derivation of the governing equations.  If we let cp vary with T and S,
!     then the governing equations will contain a lot of new terms coming from
!     gradients of cp due to spatial variation of T and S.  This complexity is
!     unnecessary due to the small variation of cp over the environmental
!     range.  In addition, the temperature T will become an implicit equation
!     of the heat, H, since H = rho(T) cp(T) T.  Note that we have also used
!     the reference density to define rho in the relation for heat:  rho(T) ->
!     rho_0.  This is known as the Boussinesq approximation.
!
!     """

   IMPLICIT NONE   
   DOUBLE PRECISION, INTENT(OUT) ::	cp
   
   cp = 3997.4

END SUBROUTINE 

SUBROUTINE compute_pressure_density(Za, Ta, Sa, Pa, comp_rho, n_dz)
!     """
!     Compute the pressure profile by integrating the density
!
!     Compute the pressure as a function of depth by integrating the density,
!     given by the temperature and salinity profiles
!	  Density profile is also returned as an array

!     Parameters
!     ----------
!     Za : ndarray
!         Array of depths in meters.
!     T : array
!         Array of temperatures (K) at the corresponding depths in `z`.
!     S : array
!         Array of salinities (psu) at the corresponding depth in `z`.
!     Returns
!     -------
!     Pa : array
!         Array of pressures (Pa) at the corresponding depth in `z`.

!	  comp_rho : array
!         Array of computed density at the corresponding depth in `z`.


	IMPLICIT NONE

	DOUBLE PRECISION, DIMENSION(n_dz), INTENT(IN) :: Za, Sa, Ta
	DOUBLE PRECISION, DIMENSION(n_dz), INTENT(INOUT) :: Pa, comp_rho	
	DOUBLE PRECISION :: g, P0, z_max, rho_sw
	INTEGER :: i, n_dz
	
	! Pressure at the sea surface
    P0 = 101325.0
    g = 9.81

	! Pressure at the surface
	Pa(1) = P0
    ! Compute the pressure at the depths in defined in z 	
    do i = 2, n_dz
		! Get the seawater density
		call density_sw(Ta(i-1), Sa(i-1), Pa(i-1), rho_sw)	
		comp_rho(i) = rho_sw
        Pa(i) = Pa(i-1) + rho_sw * g * (Za(i) - Za(i-1))
	enddo 
	! Compute the density at the maximum depth
	call density_sw(Ta(n_dz), Sa(n_dz), Pa(n_dz), rho_sw)
	comp_rho(n_dz) = rho_sw
	Write(*,*) 'End pressure Calculations'		
	
END SUBROUTINE


