program TestCABARET

	!Test program to solve the linear advection problem using CABARET
    use NetCDFModule
    implicit none
    
    INCLUDE 'netcdf.inc'
	double precision, parameter :: PI=3.14159265358979324D0

	integer, parameter :: nX=256
	integer, parameter :: nT=2000
	double precision   :: x(nX), q(nX,2)
	double precision   :: q_cell_face(nX,2)
	double precision   :: time(nT)
    double precision   :: deltaX
	double precision   :: c, CFL
	double precision   :: deltaT 
	double precision   :: std
	
	
	character(len = 8), dimension(2)	 :: dimensionNames
	integer 		   :: iX,iT
	character(len=128) :: outputFileName 
	type(netcdfFile)   :: outputNetCDFfile
	double precision   :: flux_left, flux_right
	double precision   :: q_left_extrap, q_right_extrap
	double precision   :: q_cellface_left, q_cellface_right
	double precision   :: q_half_step_left, q_half_step_centre
	
	outputFileName = 'outfile.nc'
	
	
	
	q = 1.0
	c = 1.0
	deltaX = 1.0
	CFL = 0.45
	deltaT = (deltaX * CFL)/ c
	print *, deltaT
	std  = 2.0
	x(1) = 0.0
	
	do iX=2,nX
		x(iX) = x(iX-1) + deltaX 
	enddo
	print *, 'set x'
	time(1) = 0.0
	do iT=2,nT
		time(iT) = time(iT-1) + deltaT
	enddo
	!Set initial condition
	do iX = 1,nX
	   !if (iX>nX/2) then 
		!q(iX,:) = 0.0
	   !else 
		!q(iX,:) = 1.0
	   !endif
	   
		q(iX,:) = exp(-(x(iX) - x(nX/2))*(x(iX) - x(nX/2)) / (2.0*std*std) )
	enddo 
	
	!Time Stepping Loop
	
	outputNetCDFfile = New(outputFileName)

	call WriteDimToFile(outputNetCDFfile, "x", x, "double")
	call WriteDimToFile(outputNetCDFfile, "time", time, "double")

	dimensionNames(1) =   "x"    
	dimensionNames(2) =   "time"



	call CreateNewVariable(outputNetCDFfile, 'q', dimensionNames, "double")
    
	call WriteToFile1D_Dble_SingleRecord(outputNetCDFfile, q(:,1), 'q', 1)
	
	
	!Initialise cell-face variables
	call compute_value_cell_faces(q(:,1),c,q_cell_face(:,1))

	
	do iT=2,nT
	
		
		do  iX=1,nX
			
			!Get Fluxes at cell faces x_i and x_i+1
			flux_left  = fluxes_at_cell_faces(q_cell_face(:,1),iX,c)
			flux_right = fluxes_at_cell_faces(q_cell_face(:,1),iX+1,c)
			
		
			!Take FTCS Half-Step 
			
			q_half_step_centre = q(iX,1) - ((0.5*deltaT/deltaX) * (flux_right-flux_left))
			
			
			
			!If we are on the left boundary, take a half step forward one
			!step to the left to initialise the forward extrapolation 
			!routine
			
			if (1==iX) then
				print *, 'quackulate the left half step'
				
				
				flux_left  = fluxes_at_cell_faces(q_cell_face(:,1),iX-1,c)
				flux_right = fluxes_at_cell_faces(q_cell_face(:,1), iX ,c)
				
				q_half_step_left = q(nX,1) - ((0.5*deltaT/deltaX) * (flux_right-flux_left))
				
			
			end if
			
			
			!Step the values at the cell faces forward in time by
			!linear extrapolation
			
			if (1==iX) then 
				!This step updates the cell face values at x_i
				call extrapolate_forward(q_half_step_left,  iX-1, q_cell_face)

			endif
			
			!This step updates the cell face values at x_i+1
			call extrapolate_forward(q_half_step_centre, iX,  q_cell_face)
			
			!Apply the flux limitor to the extrapolated values
			call flux_limiter_x(q_half_step_left,   iX, q_cell_face)
			call flux_limiter_x(q_half_step_centre, iX+1, q_cell_face)
			
			
			!Update the fluxes for the corrector step using the time 
			!averaged fluxes
			flux_left  = fluxes_at_cell_faces(0.5*(q_cell_face(:,2)+q_cell_face(:,1)),iX,c)
			flux_right = fluxes_at_cell_faces(0.5*(q_cell_face(:,2)+q_cell_face(:,1)),iX+1,c)
		
			!Corrector step: forward time step FTBS
			q(iX,2) = q(iX,1) - (deltaT/deltaX) * (flux_right -flux_left)
			q_half_step_left = q_half_step_centre
			
		enddo !End space loop
		
		print *, 'time step: ', iT
		call WriteToFile1D_Dble_SingleRecord(outputNetCDFfile, q(:,2), 'q', iT)
		
		q(:,1) = q(:,2)
		q(:,2) = 0.0 
		q_cell_face(:,1) = q_cell_face(:,2)
		q_cell_face(:,2) = 0.0
		
	
	enddo !END Time stepping loop
	
	
    call Delete(outputNetCDFfile)

	
contains 

subroutine compute_value_cell_faces(q,c,q_cell_face)

	double precision, dimension(:), intent(in)	:: q
	double precision, intent(in)				:: c
	double precision, dimension(:), intent(out)	:: q_cell_face
	
	integer										:: iX, nX
	
	
	nX = size(q)
	
	do iX=1,nX
		
		if (0.0>c) then
		
			q_cell_face(iX) = q(iX)
			
		else 
		
			if (1==iX) then
				q_cell_face(iX) = q(nX)
			else 
				q_cell_face(iX) = q(iX-1)
			end if
			
		end if !END 0>c
		
	
	enddo !END loop-over-x
	

end subroutine compute_value_cell_faces


!subroutine compute_value_cell_faces(q,iX,q_cellface_left,q_cellface_right)
!
!	double precision, dimension (:), intent (in) :: q
!	integer, intent(in) 		  :: iX
!	double precision, intent(out) :: q_cellface_left, q_cellface_right
!	integer 					  :: nPoints
!	
!	nPoints = size(q)
!	
!	if (1==iX) then
!		
!
!		q_cellface_left   = 0.5*(q(1) + q(nPoints))
!		q_cellface_right  = 0.5*(q(1) + q(2))
!		
!	elseif (nPoints==iX .or. 0==iX) then 
!	
!	    q_cellface_left   = 0.5*(q(nPoints) + q(nPoints-1))
!		q_cellface_right  = 0.5*(q(nPoints) + q(1))
!
!	else 
!	
!		q_cellface_left  = 0.5*(q(iX) + q(iX-1))
!		q_cellface_right = 0.5*(q(iX) + q(iX+1))
!	
!	endif
!
!	
!end subroutine compute_value_cell_faces


function fluxes_at_cell_faces(q_at_cell_face,iX,c) result(flux)
  
    double precision, dimension(:), intent(in) :: q_at_cell_face
    integer,         			    intent(in) :: iX
    double precision,				intent(in)  :: c
    
    double precision             :: flux ! output
    integer 					 :: nX
    
    nX = size(q_at_cell_face)
    
    if (iX>nX) then 
		flux = c * q_at_cell_face(1)
	
	else if (0==iX) then 
	
		flux = c * q_at_cell_face(nX)
	
	else 
		flux = c * q_at_cell_face(iX)
	endif
    
end function fluxes_at_cell_faces


subroutine flux_limiter_x(q_half_step_centre, iX, q_cell_face) 
	
	double precision, intent(in)  :: q_half_step_centre
	integer			, intent(in)  :: iX
	
	
	double precision, intent(inout) :: q_cell_face(:,:)
	
	double precision 			  :: max_limitor, min_limitor
	

	
	if (0==iX) then
		max_limitor = max(q_cell_face(nX,1), q_half_step_centre, q_cell_face(1,1)) !+ residual_flux
		min_limitor = min(q_cell_face(nX,1), q_half_step_centre, q_cell_face(1,1)) !+ residual_flux
		
		if (q_cell_face(1,2) > max_limitor) then 
			q_cell_face(1,2) = max_limitor
	 
		endif
		
		if (q_cell_face(1,2) < min_limitor) then 
			q_cell_face(1,2) = min_limitor
	
		endif
	
	
	elseif (nX<iX) then
	
		max_limitor = max(q_cell_face(nX,1), q_half_step_centre, q_cell_face(1,1)) !+ residual_flux
		min_limitor = min(q_cell_face(nX,1), q_half_step_centre, q_cell_face(1,1)) !+ residual_flux
		
		if (q_cell_face(1,2) > max_limitor) then 
			q_cell_face(1,2) = max_limitor
	 
		endif
		
		if (q_cell_face(1,2) < min_limitor) then 
			q_cell_face(1,2) = min_limitor
	
		endif
		
		
	else
		max_limitor = max(q_cell_face(iX-1,1), q_half_step_centre, q_cell_face(iX,1)) !+ residual_flux
		min_limitor = min(q_cell_face(iX-1,1), q_half_step_centre, q_cell_face(iX,1)) !+ residual_flux
	
		if (q_cell_face(iX,2) > max_limitor) then 
			q_cell_face(iX,2) = max_limitor
	 
		endif
		
		if (q_cell_face(iX,2) < min_limitor) then 
			q_cell_face(iX,2) = min_limitor
	
		endif
	
	end if
	
	
end subroutine flux_limiter_x


subroutine extrapolate_forward(q_half_step_centre, iX, q_cell_face)


	double precision, intent(in)	:: q_half_step_centre
	integer,          intent(in)	:: iX
	
	double precision, intent(inout)	:: q_cell_face(:,:)
	
	integer							:: nX
	
	nX = size(q_cell_face)
	
	if (0==iX .or. nX==iX) then 
		q_cell_face(1,2)    = (2.0 * q_half_step_centre)  - q_cell_face(nX,1)
	
	else 	
		q_cell_face(iX+1,2) = (2.0 * q_half_step_centre)  - q_cell_face(iX,1)
	
	endif 
	
end subroutine extrapolate_forward

	

end program TestCABARET
