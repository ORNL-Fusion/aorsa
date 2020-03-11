!     -----------------------------------------------------------------
!
!     global_data_m.f90
!
!
!     -----------------------------------------------------------------

MODULE global_data_m
	
	IMPLICIT NONE

! Constants

	REAL, PARAMETER :: 	pi = 3.14159265358979
	REAL, PARAMETER :: 	keV_per_erg = 1./1.6022e-9
	
	REAL, PARAMETER :: mass_e = 9.1094e-28		! mass in G
	REAL, PARAMETER :: mass_H = 1836*mass_e		! mass in G
	REAL, PARAMETER :: mass_D = 3670*mass_e		! mass in G
	REAL, PARAMETER :: mass_T = 5497*mass_e		! mass in G
	REAL, PARAMETER :: mass_3He = 5496*mass_e	! mass in G
	REAL, PARAMETER :: mass_4He = 7294*mass_e	! mass in G

! Parameters for program control

	LOGICAL :: l_read_file, l_interpolate_st_grid, l_show_graphics, l_cont_plot_f_st &
				& , l_plot_f_v, l_plot_f_cum, l_plot_f_Max, l_plot_delta_f
	
! Flags

	LOGICAL :: lHave_Distribution = .false.

! Parameters for plotting and labels

	INTEGER, PARAMETER :: kPlotParams = 6, kparam_label_length = 12, &
			& curve_label_length = 12
	REAL :: PlotParams(kPlotParams)
	CHARACTER (len=kparam_label_length) :: PlotLabels(kPlotParams), xLabel, yLabel
	CHARACTER (len=72) :: PlotTitle
	CHARACTER (len=curve_label_length) :: PlotCurveLabels(10)
	


! *******************************************************************	
	
	
END MODULE global_data_m
