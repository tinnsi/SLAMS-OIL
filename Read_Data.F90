SUBROUTINE read_measured_data()

	! This program reads data from the files and create salinity,
	! Temperature, Oxygen, density and depth data arrays
	! that can be used as input to the particle coagulation model

	USE the_info
	USE data_sim
	USE data_stations
	IMPLICIT NONE
	
	INTEGER :: AllocateStatus, i, j, ndz0, n_file
	DOUBLE PRECISION :: z0, z1, z2, para_1, para_2, para_value
	! Define array to save Density, Condictivity/Salinity, Temperature, Depth, Pressure
	DOUBLE PRECISION, DIMENSION(k_dz) :: Oa, rhoa, Sa, Ta, Za, Pa
	! Define array to save the computed Pressure, Dynamic viscosity
	DOUBLE PRECISION, DIMENSION(k_dz) :: comp_rho, dyna_vis
	
	101 FORMAT (9E15.6E2)
	1001 FORMAT(a14, 1x, a16, 1x, a16, 1x, a14, 1x, a24, 1x, a10, a14, 1x, a24, 1x, a10)
   	102 FORMAT (E15.4E2, 1X, E15.4E2, 1X, E15.4E2, 1X, E15.4E2, 1X, E15.4E2, 1X,  E15.4E2, 1X)
	1002 FORMAT(a14, 1x, a16, 1x, a16, 1x, a14, 1x, a23, 1x, a10)		
	
	write(*,*) 'Start reading measured data'
	! Open the files and create the arrays within the subroutine		
	!call data_read(file_location_data, file_name, n_skip, k_dz, n_psd, n_all, PSD, size_bin, Oa, rhoa, Sa, Ta, Za)
	call data_read(Oa, rhoa, Sa, Ta, Za)
! 	Write(*,*) 'Sa',Sa
	call compute_pressure_density(Za, Ta, Sa, Pa, comp_rho, k_dz)
	!Write(*,*) 'Sap',Sa
	! Calculate the dynamic viscosity of the seawater at Za depths
	call dynamic_viscosity_profile(Ta, Sa, Pa, dyna_vis, k_dz)
		
	
	! Get the ambient properties at different depths levels defined for the simulationand save 
	! Save thesurface conditions
	Temp_z(1) = Ta(1) - 273.15
	Sali_z(1) = Sa(1)
	Pres_z(1) = Pa(1)
	Dens_z(1) =	rhoa(1)/1000.0
	Visc_z(1) = dyna_vis(1)	* 10.0
	Depth_z(1) = Za(1)*100
	
	! Get the Energy dissipation rate in m2/s3 at 5 m depth
	CALL Energy_dissipation(Za(5), para_value)
	! Save the Energy dissipation rate in cm2/s3 
	Epsilon_z(1) = para_value*1.0e4
	! Get the average shear rate in s-1 2 m depth
	CALL Average_shear_rate(dyna_vis(1)/rhoa(1), Za(5), para_value)
	!Save the average shear rate in s-1
	Gamma_z(1) = para_value
	!Get the Kolmogorov length scale at the detph level in m
	CALL Kolmogorov_length_scale(dyna_vis(1)/rhoa(1), Epsilon_z(1)/1.0e4, para_value)
	! Save the Kolmogorov length scale in cm
	Eta_z(1) = para_value*100.0

	!Output the ambienbt parameters used for the simulations
	
!	WRITE(*,*) Sali_z(1), Temp_z(1), dens_z(1), Pres_z(1)/1e6, visc_z(1), Epsilon_z(1), &
!		 	Gamma_z(1), Depth_z(1)	
	WRITE(11,101) Sali_z(1), Temp_z(1), dens_z(1), Pres_z(1)/1e6, visc_z(1), Epsilon_z(1), &
		 	Gamma_z(1), Eta_z(1), Depth_z(1)	
	
	!write(*,*) 'ndz', ndz, 'dz', dz
	! Get the nummber of simulation depth levels in measured depth 
	ndz0 = Z_max_mea/dz	
	do i = 2, ndz0+1
		z0 = (i-1) * dz
		! Depth grid defined in the simulation
		Depth_z(i) = z0 
		j = 1
		do while (z0 .ge. Za(j-1)*100 .and. Za(j-1)*100 .lt. Z_max_mea)
			! Convert the depth to cm
			z1 = Za(j-1)*100.0
			z2 = Za(j)*100.0
			j = j + 1
! 			write(*,*) i, z1,z2, z0
		enddo
			! Get the Temperature values at dpeth levels in deg C
			para_1 = Ta(j-2) - 273.15
			para_2 = Ta(j-1) - 273.15
			CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
			Temp_z(i) =  para_value
			!write(*,*) 'para value',para_value, para_1, para_2, z0

			!Get the Salinity values at the depth levels in psu
			para_1 = Sa(j-2)
			para_2 = Sa(j-1)
			CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
			Sali_z(i) = para_value
! 			write(*,*) 'para value',para_value, para_1, para_2, z0

			!Get the Pressure values at depth levles in Pa
			para_1 = Pa(j-2)
			para_2 = Pa(j-1)
			CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
			Pres_z(i) = para_value
			!write(*,*) 'para value',para_value, para_1, para_2, z0

			!Get the Density values at depth levlesin g/cm3
			para_1 = rhoa(j-2)/1000.0
			para_2 = rhoa(j-1)/1000.0
			CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
			dens_z(i) =	para_value
			!write(*,*) 'para value',para_value, para_1, para_2, z0

			!Get the Dynamic viscosity values at depth levles in units g/cms
			para_1 = dyna_vis(j-2) * 10.0
			para_2 = dyna_vis(j-1) * 10.0
			CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
			visc_z(i) =	para_value
 			!write(*,*) 'para value',para_value, para_1, para_2, z0
			
			!Get the Energy dissipation values at depth levels in m2/s3 
			CALL Energy_dissipation(z0/100.0, para_value)
			! Save the Energy dissipation values at depth levels in cm2/s3 
			Epsilon_z(i) = para_value * 1.0e4
			
			!Get the Average shear rate values at depth levels in m2/s3 			
			CALL Average_shear_rate(visc_z(i)/10.0/(dens_z(i)*1000.0), z0/100.0, para_value)
			Gamma_z(i) = para_value
			
			!Get the Kolmogorov length scale at the detph level in m
			CALL Kolmogorov_length_scale((visc_z(i)/10.0)/(dens_z(i)*1000.0), Epsilon_z(i)/1.0e4, para_value)
			! Save the Kolmogorov length scale in cm
			Eta_z(i) = para_value*100.0
			
			!Output the ambienbt parameters used for the simulations
!			WRITE(*,*) Sali_z(i), Temp_z(i), dens_z(i), Pres_z(i)/1e6, visc_z(i), Epsilon_z(i), Gamma_z(i),Depth_z(i)
			WRITE(11,101) Sali_z(i), Temp_z(i), dens_z(i), Pres_z(i)/1.0e6, visc_z(i)/10.0,  &
				Epsilon_z(i), Gamma_z(i), Eta_z(i), Depth_z(i)
		end do

	write(*,*) 'End reading measured data'

	! Output meausre Salinity, Temperature, Density, Pressure, Dynamic Viscosity and Depth profiles
    DO i = 1, k_dz
	   !WRITE(*,*) Sa(i), Ta(i)-273.15, comp_rho(i), Pa(i)/1e6, dyna_vis(i), Za(i)
       WRITE(2,101) Sa(i), Ta(i)-273.15, comp_rho(i), Pa(i)/1e6, dyna_vis(i), Za(i)
    ENDDO
	
	! Extend the profile for deeper depth simulations
	IF (Extd_prof .eqv. .TRUE.) THEN
		CALL Extend_Profile_BM54(ndz0+1, n_file)
	ENDIF

END SUBROUTINE

SUBROUTINE data_read(O, rho, S, T, Z)
	! This subroutine reads the input data from Daly experiments and 
	! prepare an array that can be shared within the program
	
	USE data_stations
	USE data_sim	
	IMPLICIT NONE
	
	CHARACTER (len=10), DIMENSION(6) :: comment
	DOUBLE PRECISION, DIMENSION(k_dz), INTENT(INOUT) :: O, rho, S, T, Z
	DOUBLE PRECISION, DIMENSION(k_dz, n_all) :: All_Data		
	DOUBLE PRECISION, DIMENSION(k_dz) :: vol
	DOUBLE PRECISION :: mean_density, mean_salinity
	INTEGER :: i, j, k, io, AllocateStatus
	REAL, PARAMETER :: PI = 3.1415926

	All_Data(:,:) = 0.0
	vol(:) = 0.0
		
	open(UNIT=1, FILE = model_location//file_data//file_name_data,status='OLD')
		
	mean_density = 1000.0
	!write(*,*) 'read data'
		
	! Skip the header lines
	do i = 1, n_skip
		read(1,*)
	enddo

	! Get the particle size bins
	read(1,*) (comment(j),j=1,6), (agg_area(k), k=1,n_psd)
	!write(*,*) comment
 	!write(*,*) size_bin(45:53)
	!write(*,*) size_bin(1:6)

	! Skip another line 
	read(1,*)
	!write(*,*) 'start0'
	
	i = 1
	do
		!write(*,*) 'start1', comment(1), All_Data
		read(1,*, IOSTAT=io) comment(1), (All_Data(i,j), j=1,n_all)
			
		if (io /= 0) then
			EXIT
		endif
! 		write(*,*) 'start2'
		! Get the Particle numbers in #/m^3
		agg_nos(i,:) =  All_Data(i,6:5+n_psd) 
! 		write(*,*) 'agg_nos(i,:)', agg_nos(i,:)
		!Extaract the Sample volume (m3)
		vol(i) =  All_Data(i,2)
		! Extract the Oxygen (umol/kg)
		O(i) = All_Data(i,66)
		! Extract the Density (kg/m3)
		rho(i) = All_Data(i,62) + mean_density 
		! Extract the Salinity (psu/ppt)
		S(i) = All_Data(i,61)
		!write(*,*) 'S(i)',S(i)
		! Extract the temperature in degC and convert to K
		T(i) = All_Data(i,60) + 273.15
		! Extract the depth
		Z(i) = All_Data(i,1)
		!write(*,*) 'i', i
		i = i + 1	
	enddo

	do i = 1,n_psd
		! Convert the particle bin sizes assuming circular aggregates im mm
		rad_bin(i) = sqrt(agg_area(i)/pi)
	enddo

	!DEALLOCATE(All_Data)	 
	write(*,*) 'End reading data'

END SUBROUTINE

SUBROUTINE Ave_agg()
	
	USE data_stations
	USE  data_sim
	IMPLICIT NONE
	
	INTEGER :: i, j
	ave_agg_nos = 0.0
	do i= 1, n_psd
		do j= 1,nave_dz
			ave_agg_nos(i) =  ave_agg_nos(i) + agg_nos(j,i)
		enddo	
		ave_agg_nos(i) =  ave_agg_nos(i)/nave_dz
! 		write(*,*) ave_agg_nos(i)
	enddo

END SUBROUTINE

SUBROUTINE Extend_Profile_BM54(i_dx, file_no)

	!When the measured data is not available for the total depth of
	!1500m close to the DWH spill location the measured data profile
	!at BM54 values are used for the simulaitons
	!here the  new array vles are added

   USE the_info
   USE  data_sim
	IMPLICIT NONE
	CHARACTER (*), PARAMETER :: file_location_bm54 = "/Data/DWH_DATA/BM54/Output/"
	CHARACTER (*), PARAMETER :: file_name_data_tsp = "BM54_TSP.txt"
	CHARACTER (*), PARAMETER :: file_name_data_z = "BM54_Z.txt"

	DOUBLE PRECISION, DIMENSION(500, 3) :: BM54_TSP
	DOUBLE PRECISION, DIMENSION(500) :: BM54_Z
	DOUBLE PRECISION :: Z_mea_max, z0, rhosw
	DOUBLE PRECISION :: para_value, para_1, para_2, z1, z2
	INTEGER, INTENT(IN) :: i_dx
	INTEGER :: i, j, io, AllocateStatus, n_all, file_no

! 	write(*,*) file_location_bm54
	! Open the temperature, salinity and pressure data
	open(UNIT=1, FILE = model_location//file_location_bm54//file_name_data_tsp,status='OLD')

	! Skip the header lines
	read(1,*)

	! Number of columns in the file
	n_all = 3
	! Read the file
	i = 1
	do
		read(1,*, IOSTAT=io) (BM54_TSP(i,j), j=1,n_all)
		if (io /= 0) then
			EXIT
		endif
		i = i + 1
	enddo

	! Open the depth profile
	open(UNIT=2, FILE = model_location//file_location_bm54//file_name_data_z,status='OLD')
	! Skip the header lines
	read(2,*)

	! Read the file
	i = 1
	do
		read(2,*, IOSTAT=io) BM54_Z(i)
		if (io /= 0) then
			EXIT
		endif
		i = i + 1
	enddo

	101 FORMAT (8E15.6E2)
	! Get the depth at the maximum depth layer with measured data
	z0 = Depth_z(i_dx)

	! Interpolate and add the data to the existing profile
	do i = i_dx+1, ndz+1
		z0 = z0 + dz

		j = 1
		do while (z0 .ge. BM54_Z(j)*100)
			! Convert the depth to cm
			z1 = BM54_Z(j-1)*100.0
			z2 = BM54_Z(j)*100.0
			j = j + 1
 			!write(*,*) i, z1,z2, z0
		enddo

		! Update the depth array
		Depth_z(i) = z0
!		write(*,*) 'z01',z0,i
		! Get the Temperature values at dpeth levels in deg C
		para_1 = BM54_TSP(j-2,1) - 273.15
		para_2 = BM54_TSP(j-1,1) - 273.15
		CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
		Temp_z(i) =  para_value
		!write(*,*) 'para value',para_value, para_1, para_2, z0

		!Get the Salinity values at the depth levels in psu
		para_1 = BM54_TSP(j-2,2)
		para_2 = BM54_TSP(j-1,2)
		CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
		Sali_z(i) = para_value
		!write(*,*) 'para value',para_value, para_1, para_2, z0

		!Get the Pressure values at depth levles in Pa
		para_1 = BM54_TSP(j-2,3)
		para_2 = BM54_TSP(j-1,3)
		CALL interpolate(para_1,para_2,z1,z2,para_value,z0)
		Pres_z(i) = para_value

		!Calculate the corresponding density in g/cm3
		CALL density_sw(Temp_z(i)+273.15, Sali_z(i), Pres_z(i), rhosw)
		dens_z(i) =	rhosw/1000.0
		!write(*,*) 'para value',para_value, z0

		! Calculate the corresponding dynamic viscosity in units g/cms
		call dynamic_viscosity_sw(Temp_z(i)+273.15, Sali_z(i), Pres_z(i), para_value)
		! Calculate the corresponding dynamic viscosity in units kg/ms
		visc_z(i) =	para_value * 10.0
		!write(*,*) 'para value',para_value, z0

		!Get the Energy dissipation values at depth levels in m2/s3
		CALL Energy_dissipation(z0/100.0, para_value)
		! Save the Energy dissipation values at depth levels in cm2/s3
		Epsilon_z(i) = para_value * 1.0e4

		!Get the Average shear rate values at depth levels in 1/s
		CALL Average_shear_rate(visc_z(i)/10.0/(dens_z(i)*1000.0), z0/100.0, para_value)
		Gamma_z(i) = para_value

		!Get the Kolmogorov length scale at the detph level in m
		CALL Kolmogorov_length_scale((visc_z(i)/10.0)/(dens_z(i)*1000.0), Epsilon_z(i)/1.0e4, para_value)
		! Save the Kolmogorov length scale in cm
		Eta_z(i) = para_value*100.0

		!Output the ambienbt parameters used for the simulations
!		WRITE(*,*) 1, Sali_z(i), Temp_z(i), dens_z(i), Pres_z(i)/1e6, visc_z(i), Epsilon_z(i), Gamma_z(i), Depth_z(i)
		WRITE(11,101) Sali_z(i), Temp_z(i), dens_z(i), Pres_z(i)/1.0e6, visc_z(i)/10.0,  &
			Epsilon_z(i), Gamma_z(i), Depth_z(i)
	enddo

END SUBROUTINE
