#define respi 1
#define dissol 1
#define microzoo 1
#define coag2 1
#define agin 1
#define sinking 1
#define dissolved  1
#define disintegr 1
#define housekeeping 1
#define selfCollision 1
#define photo 1
#define writing	1
#define array_size 1
#define agg_comp 1
#define mass_balence 1
#define oil_added 1
#define sediment_added 1
#define Measured_Data_Input 1
#define agg_set_flux 1
#define agg_buo_flux 1
#define output_flux 1
#define mega_pic 1
#define sedimenttraps
!1
!#undef respi
!2
!#undef dissol
!3
#undef microzoo
! modify microzzop not wat when oil is there
!4
! #undef coag2 
!5
#undef agin
!6
! #undef sinking 
!7
!#undef dissolved
!8
!#undef disintegr
!9
!#undef housekeeping
!10
! #undef selfCollision
!11
#undef photo
!12
! #undef writing
!13
! #undef array_size
!14
! #undef agg_comp
!15
! #undef mass_balence
!16
! #undef oil_added
!17
! #undef sediment_added
!18
! #undef Tinna_input
!19
! #undef Measured_Data_Input
!20
!#undef agg_set_flux
!21
!#undef agg_buo_flux
!22 
#undef mega_pic
!23
!#undef sedimenttraps

MODULE the_info
	IMPLICIT NONE
   	SAVE
   
	TYPE agg_part  !Here are all the attributes of the particle (aggregate)
		!identification number, flag for aggregate at the sea floor, flag for agg or fecal pell
		! af=0: aggregate, af=1: fecal pellet
		! b=0: in watercolumn, b=1: on bottom, b=2: has been reset
      	INTEGER :: id, b, af
	  	! flag for oil only o = 1 pure oil and o = 0 coagulated particle
	  	! flag for sedi only d = 1 pure sedi and d = 0 coagulated particle
		! idl = depth level of the aggragate
		! idr = size bin of the aggaregates
      	INTEGER :: o, d, idl, idr		
		!radius, density, porosity, stickiness, depth, settl. vel
      	![um,g/cm3,-,-,cm, cm/s, #,#_total,#/n		
      	DOUBLE PRECISION :: r, rho, por, s, z, w
      	! n = number of aggregates this (super partcle/ computational particle) represents (scaling factor)
      	! p = number of primary particles in all the aggregages (ie. in one computational particle) = Nn*n
      	! Nn = number of primary particles in each (one) aggregate		
      	DOUBLE PRECISION :: n, p, Nn
		! primary particle radius um, fractal dimension, time (s)
      	DOUBLE PRECISION :: pr, frac, t
      	! minerals(calc, arag, opal, clay)		
      	DOUBLE PRECISION, DIMENSION(4) :: mineral	
      	! orgC : [molC/aggregate] orgC(:,2) - Age of the comp, orgC(1,1) - Amount of OrgC produced		  
      	DOUBLE PRECISION, DIMENSION(10,2) :: orgC, TEP
		! sedi: sediment mass [g/agg], oil; oil mass [g/agg]		
	  	DOUBLE PRECISION, DIMENSION(6,2) ::  sedi, oil
        DOUBLE PRECISION :: visc, dens  !viscosty in situ
	END TYPE
   
   
        !This is to keep track of the zooplankton
        TYPE zooplankton  
           !P and Z in molC/m3, P=phytoplankton, Z=zooplankton	   
           DOUBLE PRECISION :: P, Z
        END TYPE 

	!This is to keep track of orgC, CaCO3, etc inventory
   	DOUBLE PRECISION, DIMENSION(4) :: inventory = [0,0,0,0]
   	!OrgC, CaCO3, BSi, Clay, TEP, oil and sedi at seafloor(g)
   	DOUBLE PRECISION, DIMENSION(7) :: sf = [0,0,0,0,0,0,0] 
   	!OrgC, CaCO3, BSi, Clay, TEP, oil and sedi  net primary prod
   	DOUBLE PRECISION, DIMENSION(7) :: npp = [0,0,0,0,0,0,0] 
   	! stuff lost in wc: org,calc,opal,clay,TEP
   	DOUBLE PRECISION, DIMENSION(6) :: lost = [0,0,0,0,0,0]  
   	! mass lost in wc: org, calc1, calc2, opal, clay, TEP, oil   
   	DOUBLE PRECISION, DIMENSION(7) :: lost_mass_balnce = [0,0,0,0,0,0,0] 
   	!diameter of drops in micro meters
   	DOUBLE PRECISION, DIMENSION(6) :: drop_dia = [1.0,2.0,3.0,5.0,7.0,9.0] 
	!diameter of sedi particles in micro meters 
   	DOUBLE PRECISION, DIMENSION(6) :: sedi_dia = [1.0,2.0,3.0,5.0,7.0,9.0]	 

   	! To keep the total mass seetled at the bottom
   	DOUBLE PRECISION, DIMENSION(8) :: sf_all = [0,0,0,0,0,0,0,0]    
	
   	! Individual stickiness of the materials based on Sterling et al. (2005)
   	!stickiness of OrgC
!   	DOUBLE PRECISION :: s_OrgC = 0.5
   	DOUBLE PRECISION :: s_OrgC = 0.08 
   	!stickiness of TEP
   	DOUBLE PRECISION :: s_TEP = 0.8
   	!stickiness of oil
   	DOUBLE PRECISION :: s_oil = 0.3 
   	!stickiness of sedi 
   	DOUBLE PRECISION :: s_sedi = 0.6
   	!stickiness of bSi 
   	DOUBLE PRECISION :: s_bSi = 0.5
        !sensitivity stickiness parameter
   	DOUBLE PRECISION :: stick_param, paraT
	
   	DOUBLE PRECISION, DIMENSION(300,16) :: sedTrap
   	DOUBLE PRECISION, DIMENSION(7,7) :: ST
	
	! physical constants
	!Boltzmann's Constant J/K
   	DOUBLE PRECISION, PARAMETER :: k = 1.38d-23
   	DOUBLE PRECISION, PARAMETER :: T = 285d00 !K
   	DOUBLE PRECISION, PARAMETER :: mu = 8.9d-4
	!Shear rate s-1
   	DOUBLE PRECISION, PARAMETER :: gamma = 0.01
   	DOUBLE PRECISION, PARAMETER :: pi = 3.14159265358979323846d00
	! Gravitational acceleration (cm/s^2)	
   	DOUBLE PRECISION, PARAMETER :: g = 981 !cm/s^2
   	DOUBLE PRECISION, PARAMETER :: k_caco3_c = 5.0/86400
   	DOUBLE PRECISION, PARAMETER :: k_caco3_a = 3.2/86400
   	DOUBLE PRECISION, PARAMETER :: eta_c = 4.5
   	DOUBLE PRECISION, PARAMETER :: eta_a = 4.2
		
	! Density of seawater (g/cm3)
   	DOUBLE PRECISION, PARAMETER :: rho_sw = 1.02695d00
	! Density of Organic Carbon (g/cm3), Molar weight of Oganic Carbon (g/mol)
   	DOUBLE PRECISION, PARAMETER :: rho_orgC = 1.06d00, mw_orgC = 12.0d00!28.8d00 
   	DOUBLE PRECISION, PARAMETER :: OMCfrac = 0.54
	! Density of Calcium Carbonate (g/cm3), Molar weight of Calcium Carbonate (g/mol)	
   	DOUBLE PRECISION, PARAMETER :: rho_co = 2.8275d00, mw_co = 100.1d00	
	! Density of Biogenic silica (g/cm3), Molar weight of Biogenic silica (g/mol)
   	DOUBLE PRECISION, PARAMETER :: rho_si = 1.65d00, mw_si = 67.3d00	
	! Density of Lithogenic material (g/cm3), Molar weight of Lithogenic material (g/mol)
   	DOUBLE PRECISION, PARAMETER :: rho_li = 2.65d00, mw_li = 60d00	
	! Density of TEP (g/cm3), Molar weight of TEP (g/mol)
   	DOUBLE PRECISION, PARAMETER :: rho_TEP = 0.9, mw_TEP = 154.8d00 ! 154.8
!   	DOUBLE PRECISION, PARAMETER :: rho_TEP = 0.84, mw_TEP = 154.8d00 ! 154.8
   	!density of oil (g/cm3) 
   	DOUBLE PRECISION, PARAMETER :: rho_oil = 0.85
   	!density of sediment (g/cm3)      	
   	DOUBLE PRECISION, PARAMETER :: rho_sedi = 1.2		
		
	! Flags to keep track of the collision
   	INTEGER :: n_stick, n_collision, n_neglect
        ! Number of super-particles produced per time step
   	INTEGER :: pppts 
	
   	DOUBLE PRECISION, DIMENSION(19) :: seaBed
   	! Total depth of the water column defined in the simulation in cm
!   	INTEGER, PARAMETER :: seafloor = 63000 ! Station 1
    	INTEGER, PARAMETER :: seafloor = 150000 ! Station 2
!    	INTEGER, PARAMETER :: seafloor = 132500 ! Station 3
!    	INTEGER, PARAMETER :: seafloor = 121000 ! Station 4
!    	INTEGER, PARAMETER :: seafloor = 136000 ! Station 5
!    	INTEGER, PARAMETER :: seafloor = 159000 ! Station 6
	! Check if the available data for a shallower depth tahn the simulation depth
	! If 'yes' Extd_prof = .TRUE. and if 'no' Ext_prof = .FALSE.
! 	LOGICAL, PARAMETER :: Extd_prof = .TRUE.
   	!Thickness of a depth layer (cm)
   	INTEGER, PARAMETER :: dz = 500       
   	! Number of depth layers in the simulation
   	INTEGER, PARAMETER :: ndz = seafloor/dz 
   	! Number of total aggragates in the simulation
!    	INTEGER, PARAMETER :: a_number = 200000
   	! Depth interval that aggragate flux is calculated (cm)
   	! This should be equal to or multiplication of dz
   	INTEGER, PARAMETER :: dz_flux = 10000	!3000
   	! Time interval that aggragate flux needed to be calculated (s)
! 	INTEGER, PARAMETER :: dt_flux = 8*3600		!3*3600
   	! This should be less than the endoftime and higher than timestep 
!    	INTEGER, PARAMETER :: n_dt_flux = dt_flux/timestep
   	!Number of flux calculation depths
   	INTEGER, PARAMETER :: ndz_flux = seafloor/dz_flux
	!Number of simulation depth levels in one flux calculation depth
   	INTEGER, PARAMETER :: n_flux = dz_flux/dz
	! Depth to the middle of the intrusion layer (cm)
   	INTEGER, PARAMETER :: z_intru = 120000*0
	! Index of the depth level of the intrusion
	INTEGER, PARAMETER :: i_intr = z_intru/dz
   	! density (g/cm3), dynamic viscosity (g/cms)
   	DOUBLE PRECISION,DIMENSION(ndz+1) :: Sali_z, Temp_z, Pres_z, Dens_z, Visc_z, Depth_z
	! Average shear rate (s-1) and Energy dissipation rate (cm2/s3), Kolmogorov length scale (cm)
   	DOUBLE PRECISION,DIMENSION(ndz+1) :: Gamma_z, Epsilon_z, Eta_z			
	! The following are to view results
	! how much orgC and CaCO3 is consumed by bact and zoop at each depth 
   	DOUBLE PRECISION,DIMENSION(seafloor/dz) :: orgC_G1, orgC_b, opal_t
   	DOUBLE PRECISION,DIMENSION(seafloor/dz) :: calcC_t, calcA_t, calc_G1
   	DOUBLE PRECISION,DIMENSION(seafloor/dz) :: Th_d_234, Th_d_230    ![dpm/l]
	
END MODULE the_info	

MODULE data_sim	
	! Location of the model (user change this based on where the model is in thoer computer)
!	CHARACTER (*), PARAMETER :: model_location = "/Users/anusha/Documents/Fortran/SLAMS_TAMOC_DATA_FIX/"		
	CHARACTER (*), PARAMETER :: model_location = "/Users/tinna/git/slams/"		
	! Timestep - should not be greater than 1 day
   	INTEGER, PARAMETER :: timestep = 12*3600!*0.1	!s
   	! Number of time steps in a day
   	INTEGER, PARAMETER :: time_factor = 24*3600/timestep
	! Sumulation duration in days
	INTEGER, PARAMETER :: n_tot_days = 2000!28*6
   	! Total number of timesteps in teh simulation
   	INTEGER, PARAMETER :: endoftime = time_factor*n_tot_days*1
	! Time interval for surface loading of aggragtes of a defined size distrbution (s)
	! Choose a multiplication of the timestep
	INTEGER, PARAMETER :: dt_skip = timestep*5
	! number of time steps to skip before surface loading
   	INTEGER, PARAMETER :: n_dt_skip = dt_skip/timestep
	! Check if the available data for a shallower depth tahn the simulation depth
	! If 'yes' Extd_prof = .TRUE. and if 'no' Ext_prof = .FALSE.
	LOGICAL, PARAMETER :: Extd_prof = .TRUE.
   	! Number of total aggragates in the simulation
   	INTEGER, PARAMETER :: a_number = 200000
   	! Time interval that aggragate flux needed to be calculated (s)
	INTEGER, PARAMETER :: dt_flux = 30*24*3600   !3*3600
   	! This should be less than the endoftime and higher than timestep
   	INTEGER, PARAMETER :: n_dt_flux = dt_flux/timestep
	! Thickness of the intrusion (cm)
	INTEGER, PARAMETER :: dz_intr = 20000
	! Number of aggregates added
	INTEGER, PARAMETER :: tot_agg_intr = 100
	! (Mass flowrate can change this to later) concentration of oil into the intrusion (g/cm3)
   	DOUBLE PRECISION, PARAMETER :: intrude_oil_conc = 1.0d00/1.0d6
	! Time to initiate the intrusion (day)
	DOUBLE PRECISION, PARAMETER :: time_intrude = 4
	
	! Aggregate Size distribution
   	! Number of size bins used in aggragate radius binning
   	INTEGER, PARAMETER :: n_bins = 115
	! Array to save the size bins used in the simulations, array to save the surfave layer aggragte numbers
   	DOUBLE PRECISION, DIMENSION(n_bins) :: radi, lyr1_rad_part
	!
   	DOUBLE PRECISION, DIMENSION(2,n_bins) :: scaling	
	! To track the size, velocity, density of the particles that reach the seafloor
   	DOUBLE PRECISION, DIMENSION(n_bins) :: sizeBed, veloBed, densBed, oilBed, sediBed
	
	! to keep track of Th
	DOUBLE PRECISION :: Th_inv 
	
END MODULE data_sim



! MODULE data_stations
! 	!	Station G005
! 	! Data for the station GG1002_005S
! 	IMPLICIT NONE
! 	SAVE
!
! 	! Define the address to the folder of the data
! 	!CHARACTER (*), PARAMETER :: file_location_data = "/Users/ALD/Documents/UGA/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	CHARACTER (*), PARAMETER :: file_location_data = "/Users/anusha/Documents/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	!Define the address to the data output
! 	!CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/ALD/Documents/UGA/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/anusha/Documents/Fortran/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	! Define the filename to be opened
! 	CHARACTER (*), PARAMETER :: file_name_data = 'GG1002-005S-.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out = 'GG1002-005S-STrhoPMuZ.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out_sim = 'GG1002-005S-STrhoPMuZ-Sim.txt'
! 	! Number of the depth levels - get from the data file
! 	INTEGER, PARAMETER :: k_dz = 69 - 23
! 	! Number of particle size bins - get from the data file
! 	INTEGER, PARAMETER :: n_psd = 54
! 	! Number of header lines to be skiiped - get from the data file
! 	INTEGER, PARAMETER :: n_skip = 21
! 	! Number of total coluns to be read form the data file
! 	INTEGER, PARAMETER :: n_all = 68
! 	! Aggregate area bins for measured aggregates mm2
!   	DOUBLE PRECISION, DIMENSION(n_psd) :: agg_area
! 	! Aggregate radius bins for measured aggregates mm
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: rad_bin
! 	! Aggragate numbers in different depth layers corresponding to the sizes
! 	! in agg_area and rad_bin
! 	DOUBLE PRECISION, DIMENSION(k_dz, n_psd) :: agg_nos
! 	! Average aggragates within the defined depth
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: ave_agg_nos
! 	! No of depth layers of measured data to be avaraged
! 	INTEGER, parameter :: nave_dz = 5
! 	! Define the surface wind in m/s
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: U_wind = 5.0
! 	! Define the mixed layer dpeth in m
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: Z_mld = 50.0
! 	! Maximum Depth of the measured profile in cm
! 	! (give the rounded value to the closest simulation grid depth)
! 	INTEGER :: Z_max_mea = 4500
!
! END MODULE data_stations

! MODULE data_stations
! 	! Station G004
! 	! Data for the station GG1002_004S
! 	IMPLICIT NONE
! 	SAVE
!
! 	! Define the address to the folder of the data
! ! 	CHARACTER (*), PARAMETER :: file_location_data = "/Users/ALD/Documents/UGA/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	CHARACTER (*), PARAMETER :: file_location_data = "/Users/anusha/Documents/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	!Define the address to the data output
! ! 	CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/ALD/Documents/UGA/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/anusha/Documents/Fortran/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	! Define the filename to be opened
! 	CHARACTER (*), PARAMETER :: file_name_data = 'GG1002-004S-.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out = 'GG1002-004S-STrhoPMuZ.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out_sim = 'GG1002-004S-STrhoPMuZ-Sim.txt'
! 	! Number of the depth levels - get from the data file
! 	INTEGER, PARAMETER :: k_dz = 354 - 23
! 	! Number of particle size bins - get from the data file
! 	INTEGER, PARAMETER :: n_psd = 54
! 	! Number of header lines to be skiiped - get from the data file
! 	INTEGER, PARAMETER :: n_skip = 21
! 	! Number of total coluns to be read form the data file
! 	INTEGER, PARAMETER :: n_all = 68
! 	! Aggregate area bins for measured aggregates mm2
!   	DOUBLE PRECISION, DIMENSION(n_psd) :: agg_area
! 	! Aggregate radius bins for measured aggregates mm
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: rad_bin
! 	! Aggragate numbers in different depth layers corresponding to the sizes
! 	! in agg_area and rad_bin
! 	DOUBLE PRECISION, DIMENSION(k_dz, n_psd) :: agg_nos
! 	! Average aggragates within the defined depth
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: ave_agg_nos
! 	! No of depth layers of measured data to be avaraged
! 	INTEGER, parameter :: nave_dz = 5
! 	! Define the surface wind in m/s
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: U_wind = 5.0
! 	! Define the mixed layer dpeth in m
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: Z_mld = 50.0
! 	! Maximum Depth of the measured profile in cm
! 	! (give the rounded value to the closest simulation grid depth)
! 	INTEGER :: Z_max_mea = 33000
!
! END MODULE data_stations
!
! MODULE data_stations
! 	! Station G003
! 	! Data for the station GG1002_003S
! 	IMPLICIT NONE
! 	SAVE
! 	! Define the address to the folder of the data
! 	!CHARACTER (*), PARAMETER :: file_location_data = "/Users/ALD/Documents/UGA/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	CHARACTER (*), PARAMETER :: file_location_data = "/Users/anusha/Documents/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	!Define the address to the data output
! 	!CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/ALD/Documents/UGA/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/anusha/Documents/Fortran/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	! Define the filename to be opened
! 	CHARACTER (*), PARAMETER :: file_name_data = 'GG1002-003S-.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out = 'GG1002-003S-STrhoPMuZ.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out_sim = 'GG1002-003S-STrhoPMuZ-Sim.txt'
! 	! Number of the depth levels - get from the data file
! 	INTEGER, PARAMETER :: k_dz = 129 - 23
! 	! Number of particle size bins - get from the data file
! 	INTEGER, PARAMETER :: n_psd = 54
! 	! Number of header lines to be skiiped - get from the data file
! 	INTEGER, PARAMETER :: n_skip = 21
! 	! Number of total coluns to be read form the data file
! 	INTEGER, PARAMETER :: n_all = 68
! 	! Aggregate area bins for measured aggregates mm2
!   	DOUBLE PRECISION, DIMENSION(n_psd) :: agg_area
! 	! Aggregate radius bins for measured aggregates mm
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: rad_bin
! 	! Aggragate numbers in different depth layers corresponding to the sizes
! 	! in agg_area and rad_bin
! 	DOUBLE PRECISION, DIMENSION(k_dz, n_psd) :: agg_nos
! 	! Average aggragates within the defined depth
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: ave_agg_nos
! 	! No of depth layers of measured data to be avaraged
! 	INTEGER, parameter :: nave_dz = 5
! 	! Define the surface wind in m/s
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: U_wind = 5.0
! 	! Define the mixed layer dpeth in m
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: Z_mld = 50.0
! 	! Maximum Depth of the measured profile in cm
! 	! (give the rounded value to the closest simulation grid depth)
! 	INTEGER :: Z_max_mea = 10500
!
! END MODULE data_stations
!
! MODULE data_stations
! 	! Station G002
! 	! Data for the station GG1002_002S
! 	IMPLICIT NONE
! 	SAVE
! 	! Define the address to the folder of the data
! 	!CHARACTER (*), PARAMETER :: file_location_data = "/Users/ALD/Documents/UGA/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	CHARACTER (*), PARAMETER :: file_location_data = "/Users/anusha/Documents/Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
! 	!Define the address to the data output
! 	!CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/ALD/Documents/UGA/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	CHARACTER (*), PARAMETER :: file_location_data_output = "/Users/anusha/Documents/Fortran/SLAMS_TAMOC_DATA/Output_Input_Data/"
! 	! Define the filename to be opened
! 	CHARACTER (*), PARAMETER :: file_name_data = 'GG1002-002S-.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out = 'GG1002-002S-STrhoPMuZ.txt'
! 	CHARACTER (*), PARAMETER :: file_name_out_sim = 'GG1002-002S-STrhoPMuZ-Sim.txt'
! 	! Number of the depth levels - get from the data file
! 	INTEGER, PARAMETER :: k_dz = 174 - 23
! 	! Number of particle size bins - get from the data file
! 	INTEGER, PARAMETER :: n_psd = 54
! 	! Number of header lines to be skiiped - get from the data file
! 	INTEGER, PARAMETER :: n_skip = 21
! 	! Number of total coluns to be read form the data file
! 	INTEGER, PARAMETER :: n_all = 68
! 	! Aggregate area bins for measured aggregates mm2
!   	DOUBLE PRECISION, DIMENSION(n_psd) :: agg_area
! 	! Aggregate radius bins for measured aggregates mm
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: rad_bin
! 	! Aggragate numbers in different depth layers corresponding to the sizes
! 	! in agg_area and rad_bin
! 	DOUBLE PRECISION, DIMENSION(k_dz, n_psd) :: agg_nos
! 	! Average aggragates within the defined depth
! 	DOUBLE PRECISION, DIMENSION(n_psd) :: ave_agg_nos
! 	! No of depth layers of measured data to be avaraged
! 	INTEGER, parameter :: nave_dz = 5
! 	! Define the surface wind in m/s
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: U_wind = 7.5
! 	! Define the mixed layer dpeth in m
! 	! (assumeed for now get a value for nothern GOM in 2010)
! 	DOUBLE PRECISION, PARAMETER :: Z_mld = 50.0
! 	! Maximum Depth of the measured profile in cm
! 	! (give the rounded value to the closest simulation grid depth)
! 	INTEGER :: Z_max_mea = 15000
!
! END MODULE data_stations

MODULE data_stations
	! Station G001
	! Data for the station GG1001_001S
	IMPLICIT NONE
	SAVE

	! Location of the data (assumed the data is within the above model folder in model_location)
	CHARACTER (*), PARAMETER :: file_data = "Data/Daly_Data/GG_SD_SMP_WB_NRDA/"
	! Location of the input data profiles are saved within the model folder
	CHARACTER (*), PARAMETER :: output_data = "Output_Input_Data/"
	! Location of the model output data
	!CHARACTER (*), PARAMETER :: output_model = "Output_stn_1/"
	CHARACTER(len=16) :: output_model 
	! Define the filename to be opened
	CHARACTER (*), PARAMETER :: file_name_data = 'GG1002-001S-.txt'
	CHARACTER (*), PARAMETER :: file_name_out = 'GG1002-001S-STrhoPMuZ.txt'
	CHARACTER (*), PARAMETER :: file_name_out_sim = 'GG1002-001S-STrhoPMuZ-Sim.txt'
	! Number of the depth levels - get from the data file
	INTEGER, PARAMETER :: k_dz = 224 - 23
	! Number of particle size bins - get from the data file
	INTEGER, PARAMETER :: n_psd = 54
	! Number of header lines to be skiiped - get from the data file
	INTEGER, PARAMETER :: n_skip = 21
	! Number of total coluns to be read form the data file
	INTEGER, PARAMETER :: n_all = 68
	! Aggregate area bins for measured aggregates mm2
  	DOUBLE PRECISION, DIMENSION(n_psd) :: agg_area
	! Aggregate radius bins for measured aggregates mm
	DOUBLE PRECISION, DIMENSION(n_psd) :: rad_bin
	! Aggragate numbers in different depth layers corresponding to the sizes
	! in agg_area and rad_bin
	DOUBLE PRECISION, DIMENSION(k_dz, n_psd) :: agg_nos
	! Average aggragates within the defined depth
	DOUBLE PRECISION, DIMENSION(n_psd) :: ave_agg_nos
	! No of depth layers of measured data to be avaraged
	INTEGER, parameter :: nave_dz = 5
	! Define the surface wind in m/s
	! (assumeed for now get a value for nothern GOM in 2010)
	DOUBLE PRECISION, PARAMETER :: U_wind = 5.0
	! Define the mixed layer dpeth in m
	! (assumeed for now get a value for nothern GOM in 2010)
	DOUBLE PRECISION, PARAMETER :: Z_mld = 50.0
	! Maximum Depth of the measured profile in cm
	! (give the rounded value to the closest simulation grid depth)
	INTEGER :: Z_max_mea = 20000
END MODULE data_stations

!====================================================================

PROGRAM main

    USE the_info
    USE data_sim
    USE data_stations
    IMPLICIT NONE

    101 FORMAT (E15.4E2, 1X, E15.4E2, 1X, I5, E15.4E2, 1X, E15.4E2)
    102 FORMAT (E15.4E2, 1X, E15.4E2, 1X, I5)
    103 FORMAT (E15.4E2, 1X, I5, 1X, E15.4E2)
    104 FORMAT (E15.4E2, 1X, E15.4E2, 1X, E15.4E2)
    110 FORMAT (I4)
    111 FORMAT (E15.4E2)
    118 FORMAT (5(F12.6, 1X))
    119 FORMAT (15E18.8E2, 1X, I5)
    120 FORMAT (E15.6E2, 3X, F14.4, 1X, I5)
    125 FORMAT (E15.6E2, 3X, F14.4, 1X, I5, 1X, I5)
    121 FORMAT (E15.6E2, 1X, F14.4)
    122 FORMAT (2(E15.6E2))
!    128 FORMAT (I8, 2(E15.6E2))
    145 FORMAT (E15.4E2, 1X, E15.4E2, 1X, E15.4E2, 1X, I4)
    146 FORMAT (F14.4, 1X, 15(E15.6E2, 1X))
    147 FORMAT (5(F14.4,1X),13(E14.4E3, 1X))
    200 FORMAT (6E15.6E3, 1X, I5, 1X, I5)
    201 FORMAT (16E15.6E3, 1X, I5, 1X, I5, 1X, 3E15.6E2)
    202 FORMAT (7E15.6E3, 1X, I5)
    220 FORMAT (I5, 1X, 4E15.6E2, 1X, I5)
    240 FORMAT (5E15.6E2)
    245 FORMAT (6E15.6E2)
    250 FORMAT (I5, 1X, 4E15.6E2)
    260 FORMAT (I5, 1X, 6E15.6E2)
    270 FORMAT (2(I4), 1X, 8(E15.6E2))
    280 FORMAT (2(I4), 1X, 9(E15.6E2))
    290 FORMAT (16(E15.6E2), 1X, 2(I4))
   
    INTEGER :: ka, time, i, j, l, agg_max, z_max, temp1, temp2, tmp1, tmp2, agg_b_count
    INTEGER :: intT, intC, intD, intP, intS, intTEP, intTEP2, intCount, intO
    INTEGER :: intStick, yaff, bal1, bal2, bal3, resol, keyPulse
    INTEGER :: UNIT, free_count, loop, loop_oil, loop_sedi, loop_orgC, loop_TEP
    INTEGER :: n_dept_lyr !index for the depth layer 1 is the shallowest
    INTEGER, DIMENSION(1) :: old, seed
    INTEGER, DIMENSION(a_number), TARGET :: agg_z
    INTEGER, DIMENSION(10), TARGET :: agg_pool
    INTEGER, DIMENSION(a_number) :: free_agg
    INTEGER, POINTER :: pz
    DOUBLE PRECISION, DIMENSION(1) :: harvest
    DOUBLE PRECISION, DIMENSION(seafloor/dz,16) :: flux
    DOUBLE PRECISION, DIMENSION(seafloor/dz/5,16) :: fluxAve
    DOUBLE PRECISION :: depth, temp, dCO3c, dCO3a, max_r, organic, dummy
    DOUBLE PRECISION :: paraC, paraD, paraP, paraS, track_1, track_2
    DOUBLE PRECISION :: paraTEP, paraFrac, track_3, fudge, tf, paraOil
    DOUBLE PRECISION :: orgC1, orgC4, caco31, caco34, littleFood, bigFood, food
    DOUBLE PRECISION :: Fsize, balC, balO, balD, resolution
    DOUBLE PRECISION :: before, after, TEP_frac, year, season
    DOUBLE PRECISION :: mass, day, track4, track5, track6, track7
    DOUBLE PRECISION :: entrate_oil, entrate_sedi
    DOUBLE PRECISION, DIMENSION(seafloor/dz) :: orgCtotal, calctotal, rainratio
  	DOUBLE PRECISION, DIMENSION(seafloor/dz) :: orgCMean, calcMean, rainMean
   	DOUBLE PRECISION, DIMENSION(seafloor/dz) :: opalMean, clayMean, Th_sum_234
   	DOUBLE PRECISION, DIMENSION(seafloor/dz) :: Th_sum_230, mass_sum
   	DOUBLE PRECISION, DIMENSION(5) :: mass_i, mass_f
   	DOUBLE PRECISION, DIMENSION(15,seafloor/dz) :: data_stuff
   	!DOUBLE PRECISION, DIMENSION(n_bins,seafloor/dz) :: spec
        DOUBLE PRECISION, DIMENSION(n_bins,seafloor/dz) :: spec, spec_oil
   	DOUBLE PRECISION, DIMENSION(10,seafloor/dz) :: age
   	DOUBLE PRECISION, DIMENSION(n_bins,seafloor/dz,2,19) :: data_mega
   	DOUBLE PRECISION, DIMENSION(seafloor/dz,18) :: data_depth
   	DOUBLE PRECISION, DIMENSION(4) :: data_flux
   	DOUBLE PRECISION, DIMENSION(8) :: frac_material   
   	TYPE(agg_part), DIMENSION(a_number), TARGET :: agg
  	TYPE(agg_part), POINTER :: pa, pa1, pa2
   	TYPE(zooplankton), DIMENSION(seafloor/dz) :: mic
   	CHARACTER (len=14) :: ctime, parP, parC, parT, parS, parCount
   	CHARACTER (len=30) :: file_name, spec_name, age_name, martin_name		
	DOUBLE PRECISION :: rel_oil_mass, tot_oil_end, rel_sedi_mass, tot_sedi_end
	DOUBLE PRECISION :: rel_org_mass, tot_org_end, rel_TEP_mass, tot_TEP_end
	DOUBLE PRECISION :: rel_mineral_1, rel_mineral_2, rel_mineral_3, rel_mineral_4
	DOUBLE PRECISION :: tot_mineral_1, tot_mineral_2, tot_mineral_3, tot_mineral_4
	DOUBLE PRECISION :: min_rad
	CHARACTER*10 :: date_time_val(3)
	INTEGER :: date_time(8), i_init, r_count, nn0, n_agg_intr
	! variable to save deficit aggreagate numbers to be added to the
	DOUBLE PRECISION :: delta_n, delta_n0, delta_t, delta_s, delta_o
	!Radius of the smallest aggregate
	DOUBLE PRECISION :: init_radius
    ! Define an arrays to save the settling and buoyant flux estimates 
    ! for different components in an aggregate
     DOUBLE PRECISION,DIMENSION(ndz_flux,n_bins,8) :: set_mass_flux, buo_mass_flux


	 write(*,*) 'time step', timestep
	! Get the simulation start time and write to the screen
	call date_and_time(date_time_val(1), date_time_val(2), date_time_val(3), date_time)
	write(*,*) 'Simulation Start Date', date_time_val(1)
	write(*,*) 'Simulation Start Time', date_time_val(2)
	write(*,*) 'Simulation Start Zone', date_time_val(3)

!    READ(5,*) intT, intP, intTEP, intS, intO, intStick, intCount, keyPulse
    READ(5,*) intP, intTEP, intS, intO, intStick, intCount, keyPulse

    paraP = intP!/8.0  ! This is PP in mgC/day
    paraTEP = intTEP   ! This is TEP prod in mgC/day
    paraS = intS*1e-3  ! This is sed prod in g/day
    paraOil = intO*1e-3  ! paraOil is oil prod in g/day
    stick_param = intStick * 1e-1!1.0 / intStick 
print*, stick_param, 'stick_Param'
    WRITE(parCount,*) intCount
    parCount = adjustl(parCount)

!    WRITE(output_model,"(a13)")("INTRUSION"//trim(parCount)//"/")
    WRITE(output_model,"(a16)")("OutputSLAMS"//trim(parCount)//"/")
    print*, output_model, '  output_model', paraP, parCount
    print*, model_location, '  location'
    print*, 'Depth of seafloor' , seafloor

	intT = 20
	intC = 220
	intD = 1
	intTEP2 = 100
	
	! These are the input used by Tinna apparently for the above parameters respectively
	! 20 220 1  1000 1095 6  100 20 2  1 1 1
	write(*,*) 'time_factor', time_factor
	!Sea Surface Temperature
   	paraT = intT 
   	!CO3
   	paraC = intC  
    !zoopl. dissol. 
   	paraD = intD 	

   !!!!Add headings to all the output files
#ifdef writing

	CALL open_files()
	
#endif
	
	! Counter for the total aggregates at the bottom
	agg_b_count = 0.
	
	! Counters to check the collisions
	n_stick = 0
	n_collision = 0
	n_neglect = 0

	! Arryas to save buoyant/settlingmass fluxes at defined times
	set_mass_flux(:,:,:) = 0.0
	buo_mass_flux(:,:,:) = 0.0
	
#ifdef mass_balence
	!Initialize the parameters for mass balance 

	!Oil
	rel_oil_mass = 0.0
	tot_oil_end  = 0.0 
	
	!Sand
	rel_sedi_mass = 0.0
	tot_sedi_end = 0.0
	
	!Organic mass from primary production
	rel_org_mass = 0.0 
	tot_org_end  = 0.0
	
	!Mineral mass from primary production
	rel_mineral_1 = 0.0
	rel_mineral_2 = 0.0
	rel_mineral_3 = 0.0
	rel_mineral_4 = 0.0
	tot_mineral_1 = 0.0
	tot_mineral_2 = 0.0
	tot_mineral_3 = 0.0
	tot_mineral_4 = 0.0
	
	!TEP mass
	rel_TEP_mass = 0.0
	tot_TEP_end  = 0.0
	
#endif	
	
	! random numbers
   ka = 1
   old = 1
   seed(1) = 1834
   CALL random_seed
   CALL random_seed(size=ka)
   CALL random_seed(GET=old(1:ka))
   write(*,*) 'harvest(1)',harvest(1)
   CALL random_number(harvest(1))
   write(*,*) 'harvest(2)',harvest(1)
	
	! Get the ambient properties for the simulaiton
	CALL read_measured_data()
	
    ! Get the average particle numbers the defined depth in data_stations module
    CALL Ave_agg()	
	
!Commented out by Tinna May 8th 2019
	! Get the minimum measured aggragate radius 
	! to be considered in the model
!	r_count = 1
!	max_r = rad_bin(1)*1000.0
!	!write(*,*) 'max_r',max_r
!	DO WHILE (max_r .ge. 1.0d00)
!		max_r = max_r/1.1d00
!		r_count = r_count + 1
!	ENDDO
!	write(*,*) 'r_count', r_count
!	write(*,*) 'max_r',max_r
	
	! write the radius vector im um
   max_r = 1.0
   DO i = 1, n_bins, 1 
      radi(i) = max_r
      max_r = max_r * 1.2d00
   ENDDO 

   ! Index in the radi array 
   i_init = r_count - 1 

!	write(*,*) 'Sali_z (psu)' , Sali_z
!	write(*,*) 'Temp_z (degC)', Temp_z
!	write(*,*) 'Pres_z (MPa)', Pres_z
!	write(*,*) 'Dens_z (g/cm^3)', Dens_z
!	write(*,*) 'Dynamic Visc_z (g/cms)', Visc_z
!	write(*,*) 'Epsilon (cm2/s3)', Epsilon_z
!	write(*,*) 'Gamma (s^-1)', Gamma_z
!	write(*,*) 'Depth_z (cm)', Depth_z
! 	write(*,*) 'rad_bin (um)', rad_bin*1000
!	write(*,*) 'eta (cm)', Eta_z
   
	! Initialize variables
   time = 1
   day = 30
   !day = 365
   ! Next available spot in agg
   agg_max = 1 
   ! free_count is the counter for freed spaces in agg
   free_count = 1	
   ! This array saves the indices of the freed spaces in agg array
   free_agg = 0	
   ! An array to save the cumulative fluxes passing throuh one level
   sedTrap(:,:) = 0 
   ST(:,:) = 0 
	! This keeps track of phytoplankton					
   mic(:)%P = 0.0	
   ! This keeps track of zooplankton
   mic(:)%Z = 0.0	
	! To keep track of the materials added organic (OrgC, TEP),
	! inorganic (CaCO3,BSi, Clay)
   mass_i(:) = 0.0	
   
   seaBed = 0.0 
   lyr1_rad_part = 0.0
   year = 0

   write (*,*) 'Start of time'
   
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   
   DO WHILE (time .le. endoftime) 
! Here we do experiments where oil is not produced until day 200
     

      loop = 1
      IF (intO .eq. 0) THEN
         pppts=17
      ELSE
         pppts = 20
      ENDIF
      DO WHILE (loop .le. pppts)
!         print*, 'hello', loop, pppts
 	  !if (time .gt. 330) then
 		  !write(*,*) 'time greater than 330 first', time
	  !endif

!#ifdef Measured_Data_Input

!		IF (MOD(time-1,5) .eq. 0) THEN
      delta_n = paraP * timestep / 86400 !this is PP per timestep [mgC/ts]
      delta_s = paraS * timestep / 86400 !this is sed prod per timestep [mg/ts]
      delta_t = paraTEP * timestep / 86400 !this is TEP prod per timestep [mgC/ts]

      !Is there an oil pulse or is it constant with time? 
      ! keyPulse = 1 yes there is a pulse
      ! keyPulse = 0 no there is no pulse
      IF (keyPulse .eq. 1) THEN
!         paraOil = paraOil * 1e-1
         IF (time .le. 800*time_factor)  THEN
              delta_o = 1e-3!paraOil * timestep / 86400 * 1e-3
         ELSEIF (time .gt. 800*time_factor .and. time .le. time_factor*(887)) THEN
              delta_o = paraOil * timestep / 86400 
         ELSEIF (time .gt. time_factor*887) THEN
              delta_o = 1e-3!paraOil * timestep / 86400 * 1e-3
         ENDIF
      ELSEIF (keyPulse .eq. 0) THEN
              delta_o = paraOil * timestep / 86400
      ENDIF

!	  	  	! Primary production  aggreagates
!      	  	DO loop = 1,n_psd
!				delta_n0 = 0
!				!delta_n = part_no_GG2(loop)*1.0
!		  	  	delta_n = ave_agg_nos(loop)*nave_dz - lyr1_rad_part(loop+i_init)
!		  	  	if (time .gt. 1 .and. delta_n .gt. 1) then
!		  			delta_n = delta_n
!		  		elseif (time.eq.1 .and. delta_n .gt. 1) then
!					delta_n = ave_agg_nos(loop)*nave_dz
!  		  		elseif (time.ge.1 .and. delta_n .lt. 0) then
!  					delta_n = 0
!		  	  	endif
!
!		  		nn0 = 1
!		  	  	if (delta_n .gt. 1.0) then
!			  		delta_n0 = delta_n
!			  	  	if (delta_n .ge. nn0) then
!			  			delta_n0 = delta_n/nn0
!			  	  	else
!						delta_n0 = delta_n
!						nn0 = 1
!			  	 endif
!                                delta_n0 = delta_n0 * paraP  !production multiplication factor
!!                                print*, delta_n0, paraP, loop
!				 do i = 1,nn0
!  				
					! Assigning the aggregate index 
         			IF (free_count .gt. 1) THEN
						 ! if the agg array has a freed space use that to save the
			 			 ! particle details, free_agg array saves the indices of freed spaces in agg
            			pa => agg(free_agg(free_count-1))
            			pa%id = free_agg(free_count-1)
            			free_count = free_count-1
         			ELSEIF (free_count .eq. 1) THEN
			 			! If no freed spaces add a new location to save the particle
			 			! details (i.e.Create an empty particle)
            			pa => agg(agg_max)
            			pa%id = agg_max
            			agg_max = agg_max + 1
         			ENDIF

		 
		 			! if the total particles in the system are higher than the array defined for that
		 			! print a warning
         			IF (pa%id .gt. a_number) print*, 'END OF LIST', pa%id
		 
		 				! Create a ramdom number to be used in production subroutine
         				CALL random_number(harvest)
		
		 				! Add the aggregates
						CALL Add_Aggregate(pa, harvest, loop, delta_n, delta_t, delta_s, delta_o)

         				!CALL Data_OrgC(pa, harvest, loop, delta_n, time, paraT)
		 
		 				! Calculate the total organic matter added within the time setp. Later
		 				! this is ussed in TEP_prod to decide the TEP amount in the system
         				!organic = organic + pa%orgC(1,1)*pa%n
		 
#ifdef mass_balence		 
		 				! Keep track of the organic and mineral mass released from the primary production 
		 				rel_org_mass = rel_org_mass + pa%orgC(1,1)*pa%n * mw_orgC !g
		 				rel_mineral_1 = rel_mineral_1 + pa%mineral(1)*pa%n * mw_co	!g
		 				rel_mineral_2 = rel_mineral_2 + pa%mineral(2)*pa%n * mw_co	!g 
		 				rel_mineral_3 = rel_mineral_3 + pa%mineral(3)*pa%n * mw_si	!g
		 				rel_mineral_4 = rel_mineral_4 + pa%mineral(4)*pa%n * mw_li	!g
#endif			 


#ifdef mass_balence
		 		! Keep track of released TEP mass g
		 		rel_TEP_mass = rel_TEP_mass + pa%TEP(1,1)*pa%n*mw_TEP	
#endif


#ifdef mass_balence
		 		! Keep track of added oil mass g
		 		rel_oil_mass = rel_oil_mass + pa%oil(1,1)*pa%n
#endif


#ifdef mass_balence
		 		! Keep track of added sedi/sediment mass g
		 		rel_sedi_mass = rel_sedi_mass + pa%sedi(1,1)*pa%n
#endif
!			enddo
!		endif
!	  ENDDO
!	
!	ENDIF
	  
!#endif
         loop = loop + 1
      ENDDO !pppts
!print*, time, endoftime, npp(1)
			! Oil in the upper layers
			IF (i_intr .gt.0 .and. time .eq. time_intrude*time_factor) THEN
				n_agg_intr = 0
				! Add oil droplets in the intrusion
				DO WHILE (n_agg_intr .le. tot_agg_intr)
        			IF (free_count .gt. 1) THEN
			 			!Create an empty particle in a freed space in agg
            			pa => agg(free_agg(free_count-1))
            			pa%id = free_agg(free_count-1)
            			free_count = free_count-1
					ELSEIF (free_count .eq. 1) THEN
			 			! Create an empty particle in a new space in agg
            			pa => agg(agg_max)
            			pa%id = agg_max
            			agg_max = agg_max + 1
         			ENDIF
			 		write(*,*) 'Add oil'
					CALL random_number(harvest)
					CALL Data_Oil(pa, loop, harvest, paraT)
					n_agg_intr = n_agg_intr + 1
					
#ifdef mass_balence
! 		 		! Keep track of added oil mass g
				rel_oil_mass = rel_oil_mass + pa%oil(1,1)*pa%n
#endif
					
				ENDDO
			ENDIF


!if (time .gt. 330) then
  !write(*,*) 'time greater than 330 second', time
!endif
! End of measured data input
	  
! keeping track of how much orgC is lost to protista, zoop and bacteria
      i = 1
      track_1 = 0
      track_2 = 0
      track_3 = 0
      track4 = 0
      track5 = 0
      track6 = 0
      track7 = 0
      DO WHILE (i .le. seafloor/dz)
         track_1 = track_1 + orgC_G1(i)
         track_3 = track_3 + orgC_b(i)
         track4 = track4 + calc_G1(i)
         track6 = track6 + calcC_t(i)
         track7 = track7 + calcA_t(i)
         i = i + 1
      ENDDO

!***************************************************************************************

! find agg from same depth and put them into an array
! See if this loop around agg_max to increase the speed??	
      depth = 0
	  n_dept_lyr = 1
      DO WHILE (depth .lt. seafloor) 
         z_max = 0	! This is the maximum number of particles in the 
		 			! present depth layer 
         agg_z = 0 	! This will keep track of all the particle indices that is in the 
		 			! present depth layer 
         i = 1
         j = 1	! counter for the number of particles in the present 
		 		! depth layer 

		 ! Within this loop the z_max and agg_z for the present depth layer is caculated
         DO WHILE(i .lt. agg_max) 
            IF (agg(i)%b .eq. 0) THEN
               pz => agg_z(j)
               pa => agg(i)
               CALL check_depth(pa, pz, depth, i, j, z_max) 
            ENDIF
            i = i + 1
         ENDDO 


!i = 1
!DO WHILE (i .le. z_max)
!   IF (agg_z(i) .gt. 0) print*, depth, agg_z(i), agg(agg_z(i))%r
!   i = i + 1
!ENDDO

	! make aggregates (particles) in the same depth layer aggregate 
	! Two loops are threre to 
         i = 1		! Counter for particles
         j = 1		! Counter for particles
         DO WHILE (i .le. z_max) 
            DO WHILE (j .le. z_max) 
				! randomly find two aggregates in one depth box
  			  	CALL random_number(harvest)
				! Get the index of a random particle-1 within the layer
  				tmp1 = agg_z(FLOOR(harvest(1) * z_max))
				! If the index is zero select the partcle with the highest
				! index in the depth layer
  		  		IF (FLOOR(harvest(1)*z_max) .eq. 0) tmp1 = agg_z(z_max)

  			  	CALL random_number(harvest)
				! Get the index of a random particle-2 within the layer
  			  	tmp2 = agg_z(FLOOR(harvest(1) * z_max))
  			  	IF (FLOOR(harvest(1)*z_max) .eq. 0) tmp2 = agg_z(z_max)
				
				! Create pointers to the randomly selected particles
  			  	pa1 => agg(tmp1)
  			  	pa2 => agg(tmp2)

!This is where it is decided if the two randomly picked particles above
!do collide
				! Check if the two particles are in the water column
               IF (pa1%b .eq. 0 .and. pa2%b .eq. 0) THEN

					! Check if the two selected particles are the same 
               		IF (tmp1 .eq. tmp2) THEN
#ifdef selfColl 
                     	CALL self_collide(pa1, harvest, paraT)
#endif

                 	ELSEIF (tmp1 .ne. tmp2) THEN
#ifdef coag2
                 		CALL collide4(pa1, pa2, harvest, paraT)
#endif
                  	ENDIF
               ENDIF
               j = j + 1
            ENDDO 
!    		 	write(*,*) 'i',i
!    		 	write(*,*) 'j',j
            i = i + 1
            j = i + 1
         ENDDO 
! 		 write(*,*) 'i_final',i
! 		 write(*,*) 'j_fanal',j
		 	
! make zooplankton array using the agg_z() created above
         i = 1
         food = 0
         DO WHILE (i .le. z_max) 
            temp1 = agg_z(i)
            pa1 => agg(temp1)
            IF (pa1%af .eq. 0) THEN
               j = 1
               DO WHILE (j .le. 5) 
                  food = food + pa1%orgC(j,1)*pa1%n
                  j = j + 1
               ENDDO 
            ENDIF
            i = i + 1
         ENDDO 
		 
		 ! Find the depth cell to which the particle belong to
         CALL find_depth(depth, i)
		 ! write(*,*) 'n_dept_lyr', n_dept_lyr, i
         food = food * 1/(dz/100.0)
         mic(i)%P = food                 !molC/m3
         IF (food .gt. 0) THEN
            CALL protists(mic(i),i)
         ENDIF


#ifdef housekeeping

! Putting together the smaller particles??
         i = 1
         j = 1
         agg_pool = 0
         DO WHILE (i .le. z_max) 
            temp1 = agg_z(i)
            pa1 => agg(temp1)
            IF (pa1%Nn .eq. 1 .and. pa1%s .lt. 0.01 .and. &
                pa1%r .lt. 4) THEN
!                pa1%mineral(4) .lt. 2.8d-13) THEN
              agg_pool(j) = temp1
              j = j + 1
            ENDIF
            i = i + 1
         ENDDO 
         IF (j .gt. 3) THEN
            i = j - 2
            DO WHILE (i .ge. 1) 
               temp1 = agg_pool(j-2)
               temp2 = agg_pool(j-1)
               CALL pool(agg(temp1), agg(temp2))
               i = i - 2
               j = j - 2
            ENDDO 
         ENDIF
#endif
         depth = depth + dz
		 n_dept_lyr = n_dept_lyr + 1
      ENDDO 
	  
	  !if (time .gt. 330) then
		  !write(*,*) 'time greater than 330 third', time
	  !endif
	  
!***************************************************************************************

! Can be a separate SR like Intermed_Output
! units of sf and npp is g/m2 per 10 days, so we divide by 10 to get g/m2

         IF (MOD(time,10*time_factor) .eq. 0 .or. time.eq.endoftime) THEN
			 print*, '---------------------------', agg_max, '--------'
			 print*, 'time step:', time, 'time (hr):', time*timestep/3600.
                         print*, 'hallo Tinna', sf(1), 'temperature:', intT
			 print*, 'rel_org_mass', rel_org_mass, 'npp', npp(1)
			 print*, 'rel_oil_mass', rel_oil_mass, 'sedi:', rel_sedi_mass
			 print*, 'sf:', sf
			 print*, '-----------------*-----------------'
           	 WRITE(110,119) sf(1)/10, sf(2)/10, sf(3)/10, sf(4)/10, sf(5)/10, & 
                                sf(6)/10, sf(7)/10, &
		                npp(1)/10, npp(2)/10, npp(3)/10, npp(4)/10, npp(5)/10, &
                                npp(6)/10, npp(7)/10, stick_param, time*timestep/86400
            sf(:) = 0
            npp(:) = 0
         ENDIF
		 
!         IF (MOD(time,30*time_factor) .eq. 0 .or. time.eq.endoftime) THEN
!			 print*, '---------------------------', agg_max, '--------'
!			 print*, 'time step:', time, 'time (hr):', time*timestep/3600.
!                         print*, 'hallo Tinna', sf(1), 'temperature:', intT
!			 print*, 'rel_org_mass', rel_org_mass, 'npp', npp(1)
!			 print*, 'rel_oil_mass', rel_oil_mass, 'sedi:', rel_sedi_mass
!			 print*, 'sf:', sf
!			 print*, '-----------------*-----------------'
!           	 WRITE(110,119) sf(1)/30, sf(2)/30, sf(3)/30, sf(4)/30, sf(5)/30, & 
!                                sf(6)/30, sf(7)/30, &
!		                npp(1)/30, npp(2)/30, npp(3)/30, npp(4)/30, npp(5)/30, &
!                                npp(6)/30, npp(7)/30, stick_param, time*timestep/86400
!            sf(:) = 0
!            npp(:) = 0
!         ENDIF

#ifdef output_flux
		 
         IF (MOD(time,time_factor*30) .eq. 0) THEN
!         IF (MOD(time,time_factor*365) .eq. 0) THEN
			! calculate the flux array
			print*, 'calculating flux', time
            flux(:,:) = 0
            fluxAve(:,:) = 0
            i = 1
            DO WHILE (i .lt. agg_max)
               pa1 => agg(i)
			   ! Get the depth level l that the particle belongs to
               !CALL find_depth(pa1%z, l)
			   l = pa1%idl
			   ! Add the Organic Carbon fluxes to the array of depth levels
               j = 1
               DO WHILE (j .le. 10) 
                  flux(l,j) = flux(l,j)+pa1%orgC(j,1)*pa1%n*pa1%w/dz*timestep*time_factor
                  j = j + 1
               ENDDO 
			   ! Add the mineral fluxes to the array
               flux(l,11) = flux(l,11)+pa1%mineral(1)*pa1%n*pa1%w/dz*timestep*time_factor
               flux(l,12) = flux(l,12)+pa1%mineral(2)*pa1%n*pa1%w/dz*timestep*time_factor
               flux(l,13) = flux(l,13)+pa1%mineral(3)*pa1%n*pa1%w/dz*timestep*time_factor
               flux(l,14) = flux(l,14)+pa1%mineral(4)*pa1%n*pa1%w/dz*timestep*time_factor
			   
   				! Add the oil flux to the array
				j = 1
				DO WHILE (j .le. 6) 
				! Save the flux of oil	
   				flux(l,15) = flux(l,15)+pa1%oil(j,1)*pa1%n*pa1%w/dz*timestep*time_factor
				j = j + 1
				ENDDO	

				j = 1
				DO WHILE (j .le. 6) 
				! Save the flux of sedi	
   				flux(l,16) = flux(l,16)+pa1%sedi(j,1)*pa1%n*pa1%w/dz*timestep*time_factor
				j = j + 1
				ENDDO				
					   
               i = i + 1
            ENDDO
		
            ! average over each 50m interval or is it 40m?
            i = 1
            j = 1
            l = 1
            DO WHILE (l .le. seafloor/dz/10)
               DO WHILE (i .le. l*10)
                  DO WHILE (j .le. 16)
                     fluxAve(l,j) = fluxAve(l,j)+flux(i,j)
                     j = j + 1 
                  ENDDO
                  i = i + 1
                  j = 1
               ENDDO
               i = l*10
               l = l + 1
               j = 1
            ENDDO

            l = 1
            j = 1
            DO WHILE (l .le. seafloor/dz/10)
               DO WHILE (j .le. 16)
                  fluxAve(l,j) = fluxAve(l,j)/10.0
                  j = j + 1
               ENDDO
               j = 1
               l = l + 1
            ENDDO
			! write out the calculated flux 
   		 	! parameters for file-names
            WRITE(ctime,*) time
            ctime = adjustl(ctime)
            WRITE(parT,*) intT
            parT = adjustl(parT)
            WRITE(parC,*) intC
            parC = adjustl(parC)
            WRITE(parS,*) intS
            parS = adjustl(parS)
            WRITE(parP,*) intTEP
            parP = adjustl(parP)
			
   		! open the files to write flux into 
!       	WRITE(file_name,"(a18)")("Ave_flux_"//trim(ctime)//".dat")
!        OPEN(UNIT=102, FILE=model_location//output_model//file_name, STATUS='REPLACE', POSITION='APPEND')

		! Added by ALD 12/21/2016
!        WRITE(file_name,"(a14)")("flux_"//trim(ctime)//".dat")
!        OPEN(UNIT=5500, FILE=model_location//output_model//file_name, STATUS='REPLACE', POSITION='APPEND')			
			
!        WRITE(file_name,"(a18)")("trap_"//trim(ctime)//".dat")
!        file_name = adjustl(file_name)
!        OPEN(UNIT=53, FILE=model_location//output_model//file_name, STATUS='REPLACE', POSITION='APPEND')

#ifdef sedimenttraps
        WRITE(file_name,"(a18)")("ST_"//trim(ctime)//".dat")
        file_name = adjustl(file_name)
        OPEN(UNIT=531, FILE=model_location//output_model//file_name, STATUS='REPLACE', POSITION='APPEND')
   
   	 	! we average over a year, and return the flux in [molC/m2d]
            i = 1
            DO WHILE (i .le. seafloor/dz) 
				! Write all the data to a one file
               WRITE(550,201) flux(i,1), flux(i,2), flux(i,3),&
               flux(i,4), flux(i,5), flux(i,6), flux(i,7),&
               flux(i,8), flux(i,9), flux(i,10), flux(i,11),&
               flux(i,12), flux(i,13), flux(i,14), flux(i,15), flux(i,16), time*timestep/86400, i*dz/100, &
               (flux(i,1)+flux(i,2)+flux(i,3)+flux(i,4)+flux(i,5)+&
               flux(i,6)+flux(i,7)+flux(i,8)+flux(i,9)+flux(i,10))
			   
			   ! Write the data to different files
               WRITE(5500,201) flux(i,1), flux(i,2), flux(i,3),&
               flux(i,4), flux(i,5), flux(i,6), flux(i,7),&
               flux(i,8), flux(i,9), flux(i,10), flux(i,11),&
               flux(i,12), flux(i,13), flux(i,14), flux(i,15), flux(i,16), time*timestep/86400, i*dz/100, &
               (flux(i,1)+flux(i,2)+flux(i,3)+flux(i,4)+flux(i,5)+&
               flux(i,6)+flux(i,7)+flux(i,8)+flux(i,9)+flux(i,10))
               i = i + 1
            ENDDO
			
            flux(:,:) = 0.0
    ! write the 100m average flux 
            i = 1
            DO WHILE (i .le. seafloor/dz/10) 
               WRITE(102,201) fluxAve(i,1), fluxAve(i,2), fluxAve(i,3),&
               fluxAve(i,4), fluxAve(i,5), fluxAve(i,6), fluxAve(i,7),&
               fluxAve(i,8), fluxAve(i,9), fluxAve(i,10), fluxAve(i,11),&
               fluxAve(i,12), fluxAve(i,13), fluxAve(i,14), fluxAve(i,15), fluxAve(i,16), time*timestep/86400, i*100-50, &
               (fluxAve(i,1)+fluxAve(i,2)+fluxAve(i,3)+fluxAve(i,4)+fluxAve(i,5)+&
               fluxAve(i,6)+fluxAve(i,7)+fluxAve(i,8)+fluxAve(i,9)+fluxAve(i,10))
               i = i + 1
            ENDDO  
            fluxAve(:,:) = 0.0

   ! "measured" sediment traps
   ! we average over a year, and return the flux in [molC/m2d]
!            i = 1
!            DO WHILE (i .le. 300) 
!               WRITE(53,200) sedTrap(i,1)/day, sedTrap(i,11)/day, sedTrap(i,12)/day,&
!               sedTrap(i,14)/day, sedTrap(i,15)/day, sedTrap(i,16)/day, &
!               time*timestep/86400, i*5
!               i = i + 1
!            ENDDO
!            j = 1
!            DO WHILE (i .le. 300) 
!               WRITE(53,200) sedTrap(i,1)/day, sedTrap(i,11)/day, sedTrap(i,12)/day,&
!               sedTrap(i,14)/day, sedTrap(i,15)/day, sedTrap(i,16)/day, &
!               time*timestep/86400, j*100
!               i = i + 1
!               j = j + 1
!            ENDDO
!            sedTrap(:,:) = 0
            ST(7,:) = [20,50,100,200,500,1000,1500]
            i = 1
            DO WHILE (i .le. 7)
               WRITE(531,202) ST(1,i)/day, ST(2,i)/day, ST(3,i)/day,&
               ST(4,i)/day, ST(5,i)/day, ST(6,i)/day, ST(7,i), &
               time*timestep/86400
               i = i + 1
            ENDDO
            ST(:,:) = 0
#endif
         ENDIF
 
#endif		 
		 

!*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*^*

	! update density, velocity, ... for all particles which are in the water column
 	 i = 1 
 	!depth
 	 j = 1 
	 !sizebin
 	 l = 1
 	! Counter for the agggregates at the bottom
 	 agg_b_count = 0
 	! Number of aggregates on the surface layer
    lyr1_rad_part(:) = 0
      DO WHILE (i .lt. agg_max) 
		 ! Get all the parameters of the particles in the water column
         IF (agg(i)%b .eq. 0) THEN
 			!write(*,*) 'here-1','i, pa%z',i, pa%z
            pa => agg(i)
			! Commented by Anusha 04/207 
            !CALL find_depth(pa%z, j)
			!Get the depth level of the particle
			j = pa%idl
			 ! Commented by Anusha 04/207
			CALL find_size(pa%r, l)
			! Get the size bin of the aggregate
			l = pa%idr
			! Save the Aggregate size distribution in first depth level
			if (j .eq. 1) then
				lyr1_rad_part(l) = lyr1_rad_part(l) + pa%n
			endif
 			! Get the stickiness of the aggregate
			CALL Composite_Stickiness(pa)
			! Get the density of the aggregate
        	CALL density(pa)
			
#ifdef microzoo
            IF (mic(j)%Z .gt. 0) THEN
               CALL microzoop(pa,mic(j),harvest)
            ENDIF
#endif
#ifdef respi
			! Loos of aggregate mass due to bacterial respiration
            CALL respiration(pa, 1, dummy)
#endif 
#ifdef dissol
			! Dissolution of aggregates 
            CALL dissolution(pa)
#endif
#ifdef disintegr
			! Aggregate disintegration
            CALL disintegrate(pa, harvest)
#endif
#ifdef agin
			! Aging of aggragtes
            CALL aging(pa)
#endif
            CALL check_size(pa)

#ifdef sinking
			! Get the fractal dimension of the aggregate
			CALL fractal_dimension(pa)
			! Get the radius of the aggregate
			CALL radius(pa)
			! Get the velocity of the aggregate
			CALL velocity(pa)
!			CALL velocity_tinna_fix(pa)
!			CALL velocity_tinna(pa)
			! Settle the aggregate
			CALL sink(pa, set_mass_flux, buo_mass_flux)
			! Check if the aggregate reached the bottom
			CALL bottom(pa,l)
#endif

#ifdef dissolved
			! Check if the aggregates are completely dissolved in the water
            CALL dissolve(pa)
            CALL radius(pa)
#endif

#ifdef photo
            IF (pa%z .le. 100 .and. pa%r .gt. 1e4) THEN
               CALL photolysis(pa)
            ENDIF
#endif
			
			CALL Composite_Stickiness(pa)
			
         ENDIF
		 
 		 ! If the aggregate reached the bottom reseet the aggregate arrays
         IF (agg(i)%b .eq. 1) THEN
            free_agg(free_count) = i
            free_count = free_count + 1
			agg_b_count = agg_b_count +1
            pa => agg(i)
            CALL reset(pa)
         ENDIF
		 
		! Updatae the time of an aggragate
		CALL update_time(pa) 
		 
         i = i + 1
		 
      ENDDO 
	  
	  
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^	 

! figure out how much stuff there is at each depth
      i = 1
	  ! Saves the different types of OrgC and the minerals at diofferent levels
      data_stuff = 0
      spec = 0
      spec_oil = 0
      age = 0
      DO WHILE (i .le. seafloor/dz) 
		 ! Save the depth at each level
         data_stuff(1,i) = i*dz
         i = i + 1
      ENDDO 
	  
      i = 1
      DO WHILE (i .le. 11) 
         age(1,i) = 2**(i-1)*5*86400
         i = i + 1
      ENDDO 
	  
      i = 1
      DO WHILE (i .lt. agg_max) 
         pa => agg(i)
         IF (pa%z .lt. seafloor .and. pa%b .eq. 0) THEN
            CALL collect(pa, data_stuff, spec, spec_oil, age)
         ENDIF
         i = i + 1
      ENDDO  
	  
	  ! write how much there is into a file
      i = 1
      IF(time .gt. endoftime - 90) then
         DO WHILE (i .le. seafloor/dz) 
            orgCtotal(i) = data_stuff(2,i)+data_stuff(3,i)+data_stuff(4,i)+&
                           data_stuff(5,i)+data_stuff(6,i)+data_stuff(7,i)+&
                           data_stuff(8,i)+data_stuff(9,i)+data_stuff(10,i)+&
                           data_stuff(11,i)
            calctotal(i) = data_stuff(12,i)+data_stuff(13,i)
            orgCMean(i) = orgCMean(i) + orgCtotal(i)
            calcMean(i) = calcMean(i) + calctotal(i)
            opalMean(i) = opalMean(i) + data_stuff(14,i) 
            clayMean(i) = clayMean(i) + data_stuff(15,i)
            i = i + 1
         ENDDO 
!print*, seafloor, orgCMean(1), orgCMean(2), 'lastMonth', '*-*-*-*'
      ENDIF
	  
	  ! Calculate all the properties to make pictures
#ifdef writing
      IF (time .eq. endoftime ) THEN
         data_mega(:,:,:,:) = 0
         i = 1
         DO WHILE (i .lt. agg_max) 
            pa => agg(i)
            IF (pa%z .lt. seafloor .and. pa%b .eq. 0) THEN 
               CALL collect2(pa, data_mega)
               IF (pa%b .eq. 1) THEN 
                  print*, 'oppsie!', pa%id 
               ENDIF
            ENDIF
            i = i + 1
         ENDDO  
print*, 'ho ho'

		i = 1
        j = 1
        DO WHILE (j .le. seafloor/dz)
        	DO WHILE (i .lt. n_bins) 
				! Calculate the avrage particle density of the aggaregates
            	IF (data_mega(i,j,1,1) .gt. 0.0) THEN
                	data_mega(i,j,1,1) = data_mega(i,j,1,1)/data_mega(i,j,2,1)
               	ENDIF
			   	!Output the average density of the aggaregates
               	WRITE(20,120) data_mega(i,j,1,1), radi(i), j 
			   
			   	! Calculate the average buoyant velocity of the aggregates
               	IF (data_mega(i,j,2,2) .gt. 0.0) THEN
               		data_mega(i,j,1,2) = data_mega(i,j,1,2)/data_mega(i,j,2,2)
               	ENDIF
			   	!Output the buoyant (rising) velocities of the aggregates m/day
               	WRITE(21,120) data_mega(i,j,1,2)*864, radi(i), j 
			   	
				!Calculate the average particle settling velocity of the aggregates m/day			
               	IF (data_mega(i,j,2,15) .gt. 0.0) THEN
               		data_mega(i,j,1,15) = data_mega(i,j,1,15)/data_mega(i,j,2,15)
               	ENDIF
				!Output settling velocities of the aggregates m/day			   
               	WRITE(211,120) data_mega(i,j,1,15)*864, radi(i), j
				
				!Calculate the average stickiness of the aggregates
               	IF (data_mega(i,j,1,3) .gt. 0.0) THEN
                	data_mega(i,j,1,3) = data_mega(i,j,1,3)/data_mega(i,j,2,3)
               	ENDIF
				! Remove negative values of average stickiness (why get negative values??)
               	IF (data_mega(i,j,2,3) .le. 0.0) THEN
                	data_mega(i,j,1,3) = 0.0
               	ENDIF
			   	!Output the average stickiness of the aggregates
               	WRITE(22,120) data_mega(i,j,1,3), radi(i), j   
			   
			   	!Output the total moles of organic materials
               	data_mega(i,j,1,4) = data_mega(i,j,1,4) !/(radius(i+1)-radius(i))
               	WRITE(23,120) data_mega(i,j,1,4), radi(i), j       

			   	!Output the total moles of inorganic material
               	data_mega(i,j,1,5) = data_mega(i,j,1,5)!/(radius(i+1)-radius(i))
               	WRITE(24,120) data_mega(i,j,1,5), radi(i), j       

			   	!Output the total moles of opal
               	data_mega(i,j,1,6) = data_mega(i,j,1,6)!/(radius(i+1)-radius(i))
               	WRITE(25,120) data_mega(i,j,1,6), radi(i), j       

			   	!Output the total moles of clay
               	data_mega(i,j,1,7) = data_mega(i,j,1,7)!/(radius(i+1)-radius(i))
               	WRITE(26,120) data_mega(i,j,1,7), radi(i), j       

			   	!Calculate the averate rain ratio
               	IF (data_mega(i,j,1,8) .gt. 0.0) THEN   
                	data_mega(i,j,1,8) = data_mega(i,j,1,8)/data_mega(i,j,2,8)
               	ENDIF
				!Output the total moles of rain ratio
               	WRITE(27,120) data_mega(i,j,1,8), radi(i), j
			   
			   	!Calculate the average porosity of aggregates
               	IF (data_mega(i,j,1,9) .gt. 0.0) THEN   
                  	data_mega(i,j,1,9) = data_mega(i,j,1,9)/data_mega(i,j,2,9)
               	ENDIF
				!Output the porosity of aggregates
               	WRITE(46,120) data_mega(i,j,1,9), radi(i), j
			   
				! Output the total moles of TEP in the aggregates
               	data_mega(i,j,1,12) = data_mega(i,j,1,12)
               	WRITE(55,120) data_mega(i,j,1,12), radi(i), j
			   
			   	! Caculate the aggregate spectrum #/m^4
               	data_mega(i,j,1,16) = data_mega(i,j,1,16)/((radi(i+1)-radi(i))*1e-6)/(dz/100.0) 
			   
			   	! Output the aggragate spectrum #/m^4
               	WRITE(515,120) data_mega(i,j,1,16), radi(i), j 	
			   
			   	! Calculate the oil in each size class in each layer g/m^3
			   	data_mega(i,j,1,17) = data_mega(i,j,1,17)/dz	
			   	! Output the oil in each size class in each layer g/m^3
               	WRITE(5150,120) data_mega(i,j,1,17), radi(i), j	 
			   
			   	! oil !Number/m^4
			   	data_mega(i,j,1,18) = data_mega(i,j,1,18)/((radi(i+1)-radi(i))*1e-6)/(dz/100.0)	
               	WRITE(516,120) data_mega(i,j,1,18), radi(i), j	 			   
			   
               	i = i + 1
            ENDDO
            j = j + 1
            i = 1
         ENDDO
         !Here we are going to loop over data_mega to get a picture of how the concentration of stuff in the water
            OPEN(UNIT=160, FILE=model_location//output_model//'conc.dat', STATUS='REPLACE', POSITION='APPEND')
         !column changes with depth.
         data_depth = 0.0
         i = 1 !rad
         j = 1 !depth
         l = 1 !composition
         DO WHILE (j .le. seafloor/dz ) 
            DO WHILE (i .le. n_bins)
               DO WHILE (l .le. 18)
                  data_depth(j,l) = data_depth(j,l) + data_mega(i,j,1,l)      
                  l = l + 1
!print*, i, j, l, data_depth(j,l)
               ENDDO
               i = i + 1
               l = 1
            ENDDO
!print*, data_depth(j,1), j, i
            WRITE(160,260) j, data_depth(j,4), data_depth(j,5), data_depth(j,6), data_depth(j,12), &
                           data_depth(j,17), data_depth(j,18)
            j = j + 1
            i = 1
         ENDDO
      ENDIF
#endif


#ifdef mega_pic
   IF (MOD(time,100) .eq. 0) THEN
      print*, 'write mega_pic', spec(1,1), time
      WRITE(ctime,*) time
      ctime = adjustl(ctime)
      WRITE(spec_name,"(a18)")("spec"//trim(ctime)//".poc")
      spec_name = adjustl(spec_name)
      OPEN(UNIT=170, FILE=model_location//output_model//spec_name, STATUS='REPLACE', POSITION='APPEND')

      print*, 'write caco3_pic', spec_oil(1,1), time
      WRITE(ctime,*) time
      ctime = adjustl(ctime)
      WRITE(spec_name,"(a18)")("caco3"//trim(ctime)//".pic")
      spec_name = adjustl(spec_name)
      OPEN(UNIT=171, FILE=model_location//output_model//spec_name, STATUS='REPLACE', POSITION='APPEND')

      j = 1
      i = 1
      DO WHILE (j .le. seafloor/dz)
         DO WHILE (i .lt. n_bins)
            IF (i .eq. 1) THEN
               WRITE(170,102) spec(1,j)/radi(1), radi(i), j*dz/100
               WRITE(171,102) spec_oil(1,j)/radi(1), radi(i), j*dz/100
            ELSEIF (i .gt. 1) THEN
               WRITE(170,102) spec(i,j)/(radi(i+1)-radi(i)), radi(i), j*dz/100
               WRITE(171,102) spec_oil(i,j)/(radi(i+1)-radi(i)), radi(i), j*dz/100
            ENDIF
            i = i + 1
         ENDDO
         j = j + 1
         i = 1
      ENDDO
   ENDIF
#endif
#ifdef agg_comp
	! Output the different fractions of components of aggregates at a defined time
!      IF (MOD(time,1000) .eq. 0) THEN
!         WRITE(ctime,*) time
!         ctime = adjustl(ctime)
!         WRITE(spec_name,"(a18)")("Agg_Comp"//trim(ctime)//".dat")
!         OPEN(UNIT=172, FILE=model_location//output_model//spec_name, STATUS='REPLACE', POSITION='APPEND')
!	     i = 1
!	     DO WHILE (i .lt. agg_max)
!			j = agg(i)%idl
!			l = agg(i)%idr
!			pa => agg(i)
!			CALL Calculate_fractions(pa, frac_material)
!		  	WRITE(172,270) j, l, frac_material(1),frac_material(2),frac_material(3),frac_material(4), &
!							frac_material(5),frac_material(6),frac_material(7),frac_material(8)
!	        i = i + 1
!	     ENDDO
!	 ENDIF
#endif

#ifdef agg_set_flux
!
!	!Output the fluxes of differnt componets in aggregates at defined depths (ndz_flux)
!	! and times (n_dt_flux)
!	IF (MOD(time,n_dt_flux) .eq. 0) THEN
!   	 	WRITE(ctime,*) time*timestep/3600
!!		!write(*,*) 'time agg0', time
!   		ctime = adjustl(ctime)
!		!write(*,*) 'time agg1', time
!		file_name = adjustl("Flux_Set_")
!		! Open a file to save the settling mass fluxes
!   	 	WRITE(spec_name,*)(trim(file_name)//trim(ctime)//".dat")
!   	 	OPEN(UNIT=180, FILE=model_location//output_model//spec_name, STATUS='REPLACE', POSITION='APPEND')
!!   	 	OPEN(UNIT=180, FILE=model_location//file_location//spec_name, STATUS='REPLACE', POSITION='APPEND')
!		!write(*,*) 'time agg2', time
!		do i = 1, ndz_flux 
!			do j = 1, n_bins
!				!write(*,*) 'i,j', i,j
!				! Write the settling aggragate fluxes to a file 
!	  			WRITE(180,280) i, j, radi(j), & 
!								set_mass_flux(i,j,1), set_mass_flux(i,j,2), set_mass_flux(i,j,3), &
!								set_mass_flux(i,j,4), set_mass_flux(i,j,5), set_mass_flux(i,j,6), &
!								set_mass_flux(i,j,7), set_mass_flux(i,j,8)								
!			enddo
!		enddo
!		! Reset the array
!		set_mass_flux(:,:,:) = 0.0
!		
#endif		

#ifdef agg_buo_flux
!		WRITE(ctime,*) time*timestep/3600
!		ctime = adjustl(ctime)
!		! Open af ile to save the buoyant mass fluxes
!		file_name = adjustl("Flux_Buo_")
!   	 	WRITE(spec_name,*)(trim(file_name)//trim(ctime)//".dat")
!   	 	OPEN(UNIT=190, FILE=model_location//output_model//spec_name, STATUS='REPLACE', POSITION='APPEND')
!!   	 	OPEN(UNIT=190, FILE=file_location//spec_name, STATUS='REPLACE', POSITION='APPEND')
!		do i = 1, ndz_flux
!			do j = 1, n_bins
!				! Write the buoyant aggragate fluxes to a file
!				WRITE(190,280) i, j, radi(j), &
!							buo_mass_flux(i,j,1), buo_mass_flux(i,j,2), buo_mass_flux(i,j,3), &
!					 		buo_mass_flux(i,j,4), buo_mass_flux(i,j,5), buo_mass_flux(i,j,6), &
!					 		buo_mass_flux(i,j,7), buo_mass_flux(i,j,8)
!			enddo
!		enddo
!		! Reset the array
!		buo_mass_flux(:,:,:) = 0.0
!	ENDIF
!	
#endif

        WRITE(116,290) agg(1)%r, agg(1)%rho, agg(1)%w*864, agg(1)%s, agg(1)%Nn, agg(1)%n, agg(1)%z/100, &
                       (agg(1)%orgC(1,1)+agg(1)%orgC(2,1)+agg(1)%orgC(3,1))*mw_orgC, &
                       (agg(1)%orgC(4,1)+agg(1)%orgC(5,1)+agg(1)%orgC(6,1))*mw_orgC, &
                       (agg(1)%orgC(7,1)+agg(1)%orgC(8,1)+agg(1)%orgC(9,1)+agg(1)%orgC(10,1))*mw_orgC, &
                       (agg(1)%mineral(1)+agg(1)%mineral(2))*mw_co, &
                       agg(1)%mineral(3)*mw_si, agg(1)%sedi(1,1), agg(1)%oil(1,1), agg(1)%visc, &
                       agg(1)%dens, agg(1)%idl, agg(1)%idr
                       
        WRITE(117,290) agg(2)%r, agg(2)%rho, agg(2)%w*864, agg(2)%s, agg(2)%Nn, agg(2)%n, agg(2)%z/100, &
                       (agg(2)%orgC(1,1)+agg(2)%orgC(2,1)+agg(2)%orgC(3,1))*mw_orgC, &
                       (agg(2)%orgC(4,1)+agg(2)%orgC(5,1)+agg(2)%orgC(6,1))*mw_orgC, &
                       (agg(2)%orgC(7,1)+agg(2)%orgC(8,1)+agg(2)%orgC(9,1)+agg(2)%orgC(10,1))*mw_orgC, &
                       (agg(2)%mineral(1)+agg(2)%mineral(2))*mw_co, &
                       agg(2)%mineral(3)*mw_si, agg(2)%sedi(1,1), agg(2)%oil(1,1), agg(2)%visc, &
                       agg(2)%dens, agg(2)%idl, agg(2)%idr
                       
        WRITE(118,290) agg(3)%r, agg(3)%rho, agg(3)%w*864, agg(3)%s, agg(3)%Nn, agg(3)%n, agg(3)%z/100, &
                       (agg(3)%orgC(1,1)+agg(3)%orgC(2,1)+agg(3)%orgC(3,1))*mw_orgC, &
                       (agg(3)%orgC(4,1)+agg(3)%orgC(5,1)+agg(3)%orgC(6,1))*mw_orgC, &
                       (agg(3)%orgC(7,1)+agg(3)%orgC(8,1)+agg(3)%orgC(9,1)+agg(3)%orgC(10,1))*mw_orgC, &
                       (agg(3)%mineral(1)+agg(3)%mineral(2))*mw_co, &
                       agg(3)%mineral(3)*mw_si, agg(3)%sedi(1,1), agg(3)%oil(1,1), agg(3)%visc, & 
                       agg(3)%dens, agg(3)%idl, agg(3)%idr
                       
        WRITE(119,290) agg(4)%r, agg(4)%rho, agg(4)%w*864, agg(4)%s, agg(4)%Nn, agg(4)%n, agg(4)%z/100, &
                       (agg(4)%orgC(1,1)+agg(4)%orgC(2,1)+agg(4)%orgC(3,1))*mw_orgC, &
                       (agg(4)%orgC(4,1)+agg(4)%orgC(5,1)+agg(4)%orgC(6,1))*mw_orgC, &
                       (agg(4)%orgC(7,1)+agg(4)%orgC(8,1)+agg(4)%orgC(9,1)+agg(4)%orgC(10,1))*mw_orgC, &
                       (agg(4)%mineral(1)+agg(4)%mineral(2))*mw_co, &
                       agg(4)%mineral(3)*mw_si, agg(4)%sedi(1,1), agg(4)%oil(1,1), agg(4)%visc, & 
                       agg(4)%dens, agg(4)%idl, agg(4)%idr
                       
        !  print*, 'agg_max=', agg_max, time
        !  print*, 'agg_max=', agg_max, time
!         print*, time, agg(1)%r, agg(1)%w*864, agg(1)%z/100
      !   print*, '------------------PHYTO--------------------'
      !   print*, time, agg(13)
      !   print*, '------------------TEP--------------------'
      !   print*, time, agg(17)
      !   print*, '------------------OIL--------------------'
      !   print*, time, agg(19)
      !   print*, '------------------SEDI-------------------'
	  ! Update the time
      time = time + 1

      !we only want to see the flux for the last 400 days.
      IF (time .lt. endoftime-time_factor*200) THEN
         sizeBed(:) = 0
         veloBed(:) = 0
         densBed(:) = 0
         oilBed(:) = 0.0
         sediBed(:) = 0.0
         orgC_G1(:) = 0
         orgC_b(:) = 0
         calc_G1(:) = 0
         calcC_t(:) = 0
         calcA_t(:) = 0
         opal_t(:) = 0
         scaling(:,:) = 0
      ENDIF
   ENDDO

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
   
   i = 1
   DO WHILE (i .lt. agg_max) 
      j = 1
      DO WHILE (j .le. 10) 
         mass_f(1) = mass_f(1) + agg(i)%orgC(j,1)*agg(i)%n
         mass_f(5) = mass_f(5) + agg(i)%TEP(j,1)*agg(i)%n
         j = j + 1
      ENDDO
      mass_f(2) = mass_f(2) + (agg(i)%mineral(1)+agg(i)%mineral(2))*agg(i)%n
      mass_f(3) = mass_f(3) + agg(i)%mineral(3)*agg(i)%n
      mass_f(4) = mass_f(4) + agg(i)%mineral(4)*agg(i)%n
	  !CALL find_size(agg(i)%r, l)
	  l = agg(i)%idr
      scaling(1,l) = scaling(1,l) + 1
      scaling(2,l) = scaling(2,l) + agg(i)%n
      i = i + 1
   ENDDO
   

   WRITE(parT,*) intT
   parT = adjustl(parT)
   WRITE(parC,*) intC
   parC = adjustl(parC)
   WRITE(parP,*) intP
   parP = adjustl(parP)
   WRITE(parS,*) intTEP
   parS = adjustl(parS)
   !WRITE(file_name,"(a30)")("consu")!//trim(parT)//trim(parC)//trim(parP)
!   WRITE(file_name,*)("consu")
!   OPEN(UNIT=210, FILE=model_location//output_model//file_name, STATUS='REPLACE', POSITION='APPEND')
!   i = 1
!   DO WHILE (i .le. seafloor/dz) 
!      orgC_G1(i) = orgC_G1(i)/(dz/100)
!      orgC_b(i)  = orgC_b(i)/(dz/100)
!      calc_G1(i) = calc_G1(i)/(dz/100)
!      calcC_t(i) = calcC_t(i)/(dz/100)
!      calcA_t(i) = calcA_t(i)/(dz/100)
!      opal_t(i)  = opal_t(i)/(dz/100)
!      WRITE(210,260) i*dz/100, orgC_G1(i), orgC_b(i), &
!      calc_G1(i), calcC_t(i), calcA_t(i), opal_t(i)
!      i = i + 1
!   ENDDO

   ! Write files with the size spectrum, velocity, density etc. of stuff that falls to the seafloor.
   i = 1
   DO WHILE (i .le. n_bins)
      WRITE(111,122) radi(i), veloBed(i)
      WRITE(112,122) 0.8+0.04*(i-1), densBed(i)
      WRITE(113,122) radi(i), sizeBed(i)  !orgC [molC/m2/100days
      WRITE(114,122) radi(i), oilBed(i)
      i = i + 1
   ENDDO
   
#ifdef mass_balence

      ! Check for mass balance ALD 2/21/2017
      i = 1
      DO WHILE (i .lt. agg_max) 
         j = 1
         DO WHILE (j .le. 10) 
   		  tot_org_end = tot_org_end + agg(i)%orgC(j,1)*agg(i)%n*mw_orgC
   		  tot_TEP_end = tot_TEP_end + agg(i)%TEP(j,1)*agg(i)%n*mw_TEP
             j = j + 1
         ENDDO
	  
         j = 1
         DO WHILE (j .le. 6) 
   		  tot_oil_end = tot_oil_end + agg(i)%oil(j,1)*agg(i)%n
   		  tot_sedi_end = tot_sedi_end + agg(i)%sedi(j,1)*agg(i)%n
   		  j = j + 1
         ENDDO   
	  
   		tot_mineral_1 = tot_mineral_1 + agg(i)%mineral(1)*agg(i)%n*mw_co
   		tot_mineral_2 = tot_mineral_2 + agg(i)%mineral(2)*agg(i)%n*mw_co
   		tot_mineral_3 = tot_mineral_3 + agg(i)%mineral(3)*agg(i)%n*mw_si
   		tot_mineral_4 = tot_mineral_4 + agg(i)%mineral(4)*agg(i)%n*mw_li
      i = i + 1
      ENDDO
   
      write(*,*) 'Total Org released (g)', rel_org_mass
      write(*,*) 'OrgC in water column:', tot_org_end 
      write(*,*) 'OrgC lost to resp.', lost_mass_balnce(1)*mw_orgC
      write(*,*) 'OrgC accumul. at seafloor', sf_all(1)*mw_orgC
      write(*,*) 'Total Org in the system (g)', tot_org_end + (sf_all(1)+lost_mass_balnce(1))*mw_orgC
   
      write(*,*) 'Total TEP released (g)', rel_TEP_mass
      write(*,*) 'TEP in water column:', tot_TEP_end
      write(*,*) 'TEP lost to resp.', lost_mass_balnce(6)*mw_TEP
      write(*,*) 'TEP accumul. at seafloor', sf_all(2)*mw_TEP
      write(*,*) 'Total TEP in the system (g)', tot_TEP_end + (sf_all(2)+lost_mass_balnce(6))*mw_TEP
   
      write(*,*) 'Total Mineral 1 released (g)', rel_mineral_1
      write(*,*) 'CaCO3_c in water column:', tot_mineral_1 
      write(*,*) 'CaCO3_c lost to resp.', lost_mass_balnce(2)*mw_co
      write(*,*) 'CaCO3_c accumul. at seafloor', sf_all(3)*mw_co
      write(*,*) 'Total Mineral 1 in the system (g)', tot_mineral_1 + (sf_all(3)+lost_mass_balnce(2))*mw_co
   
      write(*,*) 'Total Mineral 2 released (g)', rel_mineral_2
      write(*,*) 'CaCO3_a in water column:', tot_mineral_2
      write(*,*) 'CaCO3_a lost to resp.', lost_mass_balnce(3)*mw_co
      write(*,*) 'CaCO3_a accumul. at seafloor', sf_all(4)*mw_co
      write(*,*) 'Total Mineral 2 in the system (g)', tot_mineral_2 + (sf_all(4)+lost_mass_balnce(3))*mw_co
   
      write(*,*) 'Total Mineral 3 released (g)', rel_mineral_3
      write(*,*) 'Opal in water column:', tot_mineral_3
      write(*,*) 'Opal lost to resp.', lost_mass_balnce(4)*mw_si
      write(*,*) 'Opal accumul. at seafloor', sf_all(5)*mw_si
      write(*,*) 'Total Mineral 3 in the system (g)', tot_mineral_3 + (sf_all(5)+lost_mass_balnce(4))*mw_si   

      write(*,*) 'Total Mineral 4 released (g)', rel_mineral_4   
      write(*,*) 'Total Mineral 4 in the system (g)', tot_mineral_4 + (sf_all(6)+lost_mass_balnce(5))*mw_li   

      write(*,*) 'Total Oil released (g)', rel_oil_mass   
      write(*,*) 'Oil in water column:', tot_oil_end 
      write(*,*) 'Oil lost to resp.', lost_mass_balnce(7)
      write(*,*) 'Oil accumul. at seafloor', sf_all(7)
      write(*,*) 'Total Oil in the system (g)', tot_oil_end + sf_all(7) + lost_mass_balnce(7)

      write(*,*) 'Total Sand released (g)', rel_sedi_mass   
      write(*,*) 'Sediments accumul. at seafloor', sf_all(8)
      write(*,*) 'Total Sand in the system (g)', tot_sedi_end + sf_all(8)  

#endif   
   
   write(*,*) 'n_stick', n_stick
   write(*,*) 'n_collision', n_collision
   write(*,*) 'percentage of successful collision within stick2 from total' , n_stick*100/n_collision
   
   write(*,*)  'n_neglect', n_neglect
  
#ifdef respi
print*, 'Respiration ON'
#endif
#ifdef dissol
print*, 'Dissolution ON'
#endif
#ifdef microzoo
print*, 'Microzooplankton ON'
#endif
#ifdef coag2
print*, 'Coagulation ON'
#endif
#ifdef agin
print*, 'Aging ON'
#endif
#ifdef sinking
print*, 'Sinking ON'
#endif
#ifdef dissolved  
print*, 'Dissolve ON'
#endif
#ifdef disintegr 
print*, 'Disintegrate ON'
#endif
#ifdef housekeeping 
print*, 'Housekeeping ON'
#endif
#ifdef selfCollision 
print*, 'Self collision ON'
#endif
#ifdef photo 
print*, 'Photo ON'
#endif
#ifdef writing	
print*, 'Writing ON'
#endif
#ifdef array_size 
print*, 'Array_size ON'
#endif
#ifdef agg_comp 
print*, 'Agg_comp ON'
#endif
#ifdef mass_balence 
print*, 'Mass balance ON'
#endif
#ifdef oil_added 
print*, 'Oil added ON'
#endif
#ifdef sediment_added 
print*, 'Sediment added ON'
#endif
#ifdef Measured_Data_Input 
print*, 'Measured Data Input ON'
#endif



	! Get the simulation End time and write to the screen
	call date_and_time(date_time_val(1), date_time_val(2), date_time_val(3), date_time)
	write(*,*) 'Simulation End Date', date_time_val(1)
	write(*,*) 'Simulation End Time', date_time_val(2)
	write(*,*) 'Simulation End Zone', date_time_val(3)

END PROGRAM
