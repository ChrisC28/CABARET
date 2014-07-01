program TestCABARET2D

	!Test program to solve the linear advection problem using CABARET
    use NetCDFModule
    implicit none
    
    INCLUDE 'netcdf.inc'
	double precision, parameter :: PI=3.14159265358979324D0

	integer, parameter :: nX=128,nY=128
	integer, parameter :: nT=500
	double precision   :: x(nX), y(nY), q(nX,nY,2)
	double precision   :: q_cellface_x(nX,nY,2), q_cellface_y(nX,nY,2) 
	double precision   :: c(2)
	double precision   :: time(nT)
    double precision   :: deltaX, deltaY
	double precision   ::  CFL
	double precision   :: deltaT 
	double precision   :: std
	
	
	character(len = 8), dimension(3)	 :: dimensionNames
	integer 		   :: iX,iY, iT
	character(len=128) :: outputFileName 
	type(netcdfFile)   :: outputNetCDFfile
	double precision   :: flux_left, flux_right
	double precision   :: flux_south, flux_north
	
	double precision   :: q_left_extrap,    q_right_extrap
	double precision   :: q_north_extrap,   q_south_extrap
	
	double precision   :: q_cellface_left,  q_cellface_right
	double precision   :: q_cellface_south, q_cellface_north

	double precision   :: q0, q_half_step_left
	double precision   :: q_half_step_centre, q_half_step_south(nX)
	double precision   :: oneOnDeltaX, oneOnDeltaY
	
	outputFileName = 'linearAdvection2D.nc'
	 
	
	
	c(1) = 1.0
	c(2) = 1.0
	deltaX = 1.0
	deltaY = 1.0
	
	
	oneOnDeltaX = 1.0/deltaX
	oneOnDeltaY = 1.0/deltaY
	
	
	CFL = 0.4
	deltaT = (deltaX * CFL)/ c(1)
	print *, deltaT
	std  = 4.0
	x(1) = 0.0
	y(1) = 0.0
	do iX=2,nX
		x(iX) = x(iX-1) + deltaX 
	enddo
	do iY=2,nY
		y(iY) = y(iY-1) + deltaY 
	enddo
	
	
	time(1) = 0.0
	do iT=2,nT
		time(iT) = time(iT-1) + deltaT
	enddo
	
	!Set initial condition
	do iY = 1,nY
		q0 = exp(-(y(iY) - y(nY/2))**2 / (2.0*std*std) )
	
		do iX = 1,nX
	   
			q(iX,iY,1) = q0*exp(-(x(iX) - x(nX/2))**2 / (2.0*std*std) ) 
	
		enddo 
	enddo
	
	!Time Stepping Loop
	
	outputNetCDFfile = New(outputFileName)

	call WriteDimToFile(outputNetCDFfile, "x", x, "double")
	call WriteDimToFile(outputNetCDFfile, "y", y, "double")

	call WriteDimToFile(outputNetCDFfile, "time", time, "double")

	dimensionNames(1) =   "x"
	dimensionNames(2) =   "y"    
    dimensionNames(3) =   "time"
	


	call CreateNewVariable(outputNetCDFfile, 'q', dimensionNames, "double")
    print *, 'Writing initial output'
	call WriteToFile2D_Dble_SingleRecord(outputNetCDFfile, q(:,:,1), 'q', 1)
	print *, 'output written!'
    !Initialise the variables at cell boundaries
    call compute_value_cell_faces(q(:,:,1),c,q_cellface_x(:,:,1),q_cellface_y(:,:,1))
    
    
	
	do iT=2,nT
	
	 do iY=1,nY
	  do iX=1,nX
	    !==========================================================!
	    !FTCS half-step
	    !==========================================================!
			
	    !Compute fluxes in the x direction
	    
	    !Guess values of the state variable the cell faces to the left 
	    !and right 
	   

	    !Use the values of the state variables to compute the fluxes 
	    !through the left and right cell faces
	   
		flux_left  = fluxes_at_cell_faces(q_cellface_x(:,iY,1), iX,   c(1))
	    flux_right = fluxes_at_cell_faces(q_cellface_x(:,iY,1), iX+1, c(1))
		
		!Use the values of the state variables to compute the fluxes 
	    !through the south and north cell faces
		flux_south =  fluxes_at_cell_faces(q_cellface_y(iX,:,1), iY,  c(2))
		flux_north =  fluxes_at_cell_faces(q_cellface_y(iX,:,1), iY+1,c(2))
	   
	    !Take the half step at the centered point
		q_half_step_centre = q(iX,iY,1) - 0.5*deltaT*(              & 
		                     oneOnDeltaX * (flux_right-flux_left) + &
		                     oneOnDeltaY * (flux_north-flux_south)  )
		

		if (1==iX) then
		
			!If we are on the left boundary, we need to compute the forward
			!time half step to the left, using the periodic bc's
			!in order to extrapolate the state variable forward in time
			
			!Guess values of the state variable the cell faces to the left 
			!and right 
	    
			!Use the values of the state variables to compute the fluxes 
			!through the left and right cell faces
			
		
			flux_left  = fluxes_at_cell_faces(q_cellface_x(:,iY,1), iX-1, c(1))
			flux_right = fluxes_at_cell_faces(q_cellface_x(:,iY,1), iX,   c(1))
		
			!Use the values of the state variables to compute the fluxes 
			!through the south and north cell faces
			flux_south =  fluxes_at_cell_faces(q_cellface_y(nX,:,1), iY,     c(2))
			flux_north =  fluxes_at_cell_faces(q_cellface_y(nX,:,1), iY+1   ,c(2))
			
			!Take the half step at the point to the left

			q_half_step_left = q(nX,iY,1) - 0.5*deltaT*(                & 
							   oneOnDeltaX * (flux_right-flux_left) +   &
							   oneOnDeltaY * (flux_north-flux_south)    )
			 
		end if
		
		if (1==iY) then
		
			!If we are on the southern boundary, we need to compute the forward
			!time half step to the north, using the periodic bc's
			!in order to extrapolate the state variable forward in time
			
			!Guess values of the state variable the cell faces to the left 
			!and right 
	    
			flux_left  = fluxes_at_cell_faces(q_cellface_x(:,nY,1), iX,   c(1))
			flux_right = fluxes_at_cell_faces(q_cellface_x(:,nY,1), iX+1, c(1))
		
			!Use the values of the state variables to compute the fluxes 
			!through the south and north cell faces
			flux_south =  fluxes_at_cell_faces(q_cellface_y(iX,:,1), iY-1, c(2))
			flux_north =  fluxes_at_cell_faces(q_cellface_y(iX,:,1), iY   ,c(2))
			
			
			q_half_step_south(iX) = q(iX,nY,1) - 0.5*deltaT*(            & 
								    oneOnDeltaX * (flux_right-flux_left) + &
								    oneOnDeltaY * (flux_north-flux_south) )
		end if
		
		
		!Extrapolate from t=n+1/2 to t=n+1 in the x direction
	    !call compute_value_cell_faces(q(:,iY,1), iX-1, q_cellface_left,q_cellface_right)
		
		!Step the values at the cell faces forward in time by
		!linear extrapolation
			
		if (1==iX) then 
		!	!This step steps forward in time the cell face values at y_i,x_i
			call extrapolate_forward(q_half_step_left,  iX-1, q_cellface_x(:,iY,:))
			call flux_limitor_x2(q_cellface_x, q_cellface_y, q_half_step_left, iX, iY, deltaX, deltaY, deltaT, c(2))
		endif
		
		if (1==iY) then
		!	!This step steps forward in time the cell face values at y_i,x_i
			call extrapolate_forward(q_half_step_south(iX), iY-1, q_cellface_y(iX,:,:))
			call flux_limitor_y2(q_cellface_y, q_cellface_x, q_half_step_south(iX), iX, iY, deltaX, deltaY, deltaT, c(1))
		endif	
		
		!This step updates the cell face values at x_i+1
		
		!call flux_limitor_x2(q_cellface_x, q_cellface_y, q_half_step_left, iX-1, iY, deltaX, deltaY, deltaT, c(2))
		call extrapolate_forward(q_half_step_centre, iX,  q_cellface_x(:,iY,:))
		call extrapolate_forward(q_half_step_centre, iY,  q_cellface_y(iX,:,:))
		
		call flux_limitor_x2(q_cellface_x, q_cellface_y, q_half_step_centre, iX+1, iY, deltaX, deltaY, deltaT, c(2))
		call flux_limitor_y2(q_cellface_y, q_cellface_x, q_half_step_centre, iX, iY+1, deltaX, deltaY, deltaT, c(1))

		
		
		!
		
		flux_left  = fluxes_at_cell_faces(0.5*(q_cellface_x(:,iY,1)+q_cellface_x(:,iY,2)), iX,   c(1))
	    flux_right = fluxes_at_cell_faces(0.5*(q_cellface_x(:,iY,1)+q_cellface_x(:,iY,2)), iX+1, c(1))
		
		!Use the values of the state variables to compute the fluxes 
	    !through the south and north cell faces
		flux_south =  fluxes_at_cell_faces(0.5*(q_cellface_y(iX,:,1)+q_cellface_y(iX,:,2)), iY,  c(2))
		flux_north =  fluxes_at_cell_faces(0.5*(q_cellface_y(iX,:,1)+q_cellface_y(iX,:,2)), iY+1,c(2))
			
		
		!Finally, take another BTCS step
		q(iX,iY,2) = q(iX,iY,1)   - deltaT*( & 
		              oneOnDeltaX * (flux_right-flux_left) + &
		              oneOnDeltaY * (flux_north-flux_south) )
		
		 q_half_step_left  = q_half_step_centre
		 
	  enddo !END iX
		q_half_step_south = q_half_step_centre
		
	 enddo !END iY
		
	
	  
	  print *, 'time step: ', iT
	  call WriteToFile2D_Dble_SingleRecord(outputNetCDFfile, q(:,:,2), 'q', iT)

	  q(:,:,1)  = q(:,:,2)
	  
	  q_cellface_x(:,:,1) = q_cellface_x(:,:,2)
	  q_cellface_x(:,:,2) = 0.0
	  q_cellface_y(:,:,1) = q_cellface_y(:,:,2)
	  q_cellface_y(:,:,2) = 0.0
	enddo !!END Time Stepping routine
	
	call Delete(outputNetCDFfile)



contains

subroutine compute_value_cell_faces(q,c,q_cellface_x,q_cellface_y)

	double precision, intent(in)	:: q(:,:)
	double precision, intent(in)	:: c(2)
	
	double precision, dimension (:,:), intent(out)	:: q_cellface_x,q_cellface_y
	
	integer 						:: iX, iY, nX, nY
	
	nX = size(q,1)
	nY = size(q,2)
	
	do iY=1,nY 
		do iX=1,nX
		
			if (0.0>c(1)) then
				
				q_cellface_x(iX,iY) = q(iX,iY)
				
			else
				
				if (1==iX) then 
					q_cellface_x(iX,iY) = q(nX,iY)
				else
					q_cellface_x(iX,iY) = q(iX-1,iY)
				endif !END 1==iX
				
			endif !END 0>c(1)
			
			if (0.0>c(2)) then
				
				q_cellface_y(iX,iY) = q(iX,iY)

			else	
				
				if (1==iY) then
					q_cellface_y(iX,iY) = q(iX,nY)
				else 
					q_cellface_y(iX,iY) = q(iX,iY-1)	
				end if !END 1==iY
				
				
			end if !0>c(2)
			
	
		enddo !END x-loop
	enddo  !END y-loop


end subroutine compute_value_cell_faces

function fluxes_at_cell_faces(q_at_cell_face, iX, c) result(flux)
  
    double precision, dimension(:), intent(in) :: q_at_cell_face
    integer,         			    intent(in) :: iX
    double precision,				intent(in) :: c
    
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


subroutine flux_limitor_x2(q_cellface_x, q_cellface_y, q_half_step, iX, iY, deltaX, deltaY,deltaT, v)

	double precision, intent(inout) 	:: q_cellface_x(:,:,:)
	
	double precision, intent(in)		:: q_cellface_y(:,:,:)
	double precision, intent(in)		:: q_half_step
	integer, intent(in)					:: iX, iY
	double precision, intent(in)		:: deltaX, deltaY, deltaT, v
	
	double precision					:: max_limiter
	double precision					:: min_limiter
	integer								:: nX, nY
	double precision					:: source_term
	
	nX = size(q_cellface_x,1)
	nY = size(q_cellface_x,2)
	
	
	if (iY == 0) then 
		source_term = -1.0 * (v / deltaY) * ( q_cellface_y(iX,nY,1)  - q_cellface_y(iX,nY-1,1) ) * deltaT
	elseif (iY==1) then
		source_term = -1.0 * (v / deltaY) * ( q_cellface_y(iX,1,1)  - q_cellface_y(iX,nY,1) ) * deltaT
	else 
		source_term = -1.0 * (v / deltaY) * ( q_cellface_y(iX,iY,1) - q_cellface_y(iX,iY-1,1) ) * deltaT
	endif
	
	if (0==iX) then
		
		max_limiter = max(q_cellface_x(nX-1,iY,1),q_cellface_x(nX,iY,1), q_half_step) + source_term
		min_limiter = min(q_cellface_x(nX-1,iY,1),q_cellface_x(nX,iY,1), q_half_step) + source_term
		
		if (q_cellface_x(nX,iY,2) > max_limiter) then
			
			q_cellface_x(nX,iY,2) = max_limiter
		
		endif
		
		if (q_cellface_x(nX,iY,2) < min_limiter) then
			
			q_cellface_x(nX,iY,2) = min_limiter
		
		endif
		
	
	else if (1==iX) then 
		
		max_limiter = max(q_cellface_x(nX,iY,1),q_cellface_x(1,iY,1),q_half_step) + source_term
		min_limiter = min(q_cellface_x(nX,iY,1),q_cellface_x(1,iY,1),q_half_step) + source_term
		
		
		if (q_cellface_x(1,iY,2) > max_limiter) then
			
			q_cellface_x(1,iY,2) = max_limiter
		
		endif
		
		if (q_cellface_x(1,iY,2) < min_limiter) then
			
			q_cellface_x(1,iY,2) = min_limiter
		
		endif
	
	else 
	
		max_limiter = max(q_cellface_x(iX-1,iY,1), q_cellface_x(iX,iY,1), q_half_step) + source_term
		min_limiter = min(q_cellface_x(iX-1,iY,1), q_cellface_x(iX,iY,1), q_half_step) + source_term
		
		if (q_cellface_x(iX,iY,2) > max_limiter) then
			
			q_cellface_x(iX,iY,2) = max_limiter
		
		endif
		
		if (q_cellface_x(iX,iY,2) < min_limiter) then
			
			q_cellface_x(iX,iY,2) = min_limiter
		
		endif	
	endif


end subroutine flux_limitor_x2


subroutine flux_limitor_y2(q_cellface_y, q_cellface_x, q_half_step, iX, iY, deltaX, deltaY,deltaT, u)

	double precision, intent(inout) 	:: q_cellface_y(:,:,:)
	
	double precision, intent(in)		:: q_cellface_x(:,:,:)
	double precision, intent(in)		:: q_half_step
	integer, intent(in)					:: iX, iY
	double precision, intent(in)		:: deltaX, deltaY, deltaT, u
	
	double precision					:: max_limiter
	double precision					:: min_limiter
	integer								:: nX, nY
	double precision					:: source_term
	
	nX = size(q_cellface_x,1)
	nY = size(q_cellface_x,2)
	
	
	if (iX == 0) then 
		source_term = - 2.0 * (u / deltaX) * (q_cellface_x(nX,iY,1) -  q_cellface_x(nX-1,iY,1)) * deltaT
	elseif (iX==1) then
		source_term = - 2.0 * (u / deltaX) * (q_cellface_x(1,iY,1)  -  q_cellface_x(nX,iY,1)) * deltaT
	else 
		source_term = - 2.0 * (u / deltaX) * (q_cellface_x(iX,iY,1) -  q_cellface_x(iX-1,iY,1)) * deltaT
	endif
	
	if (0==iY) then
		
		max_limiter = max(q_cellface_y(iX,nY-1,1),q_cellface_y(iX,nY,1), q_half_step) + source_term
		min_limiter = min(q_cellface_y(iX,nY-1,1),q_cellface_y(iX,nY,1), q_half_step) + source_term
		
		if (q_cellface_y(iX,nY,2) > max_limiter) then
			
			q_cellface_y(iX,nY,2) = max_limiter
		
		endif
		
		if (q_cellface_y(iX,nY,2) < min_limiter) then
			
			q_cellface_y(iX,nY,2) = min_limiter
		
		endif
		
	
	else if (1==iY) then 
		
		max_limiter = max(q_cellface_y(iX,1,1),q_cellface_x(iX,1,1),q_half_step) + source_term
		min_limiter = min(q_cellface_y(iX,1,1),q_cellface_x(iX,1,1),q_half_step) + source_term
		
		
		if (q_cellface_y(iX,1,2) > max_limiter) then
			
			q_cellface_y(iX,1,2) = max_limiter
		
		endif
		
		if (q_cellface_y(iX,1,2) < min_limiter) then
			
			q_cellface_y(iX,1,2) = min_limiter
		
		endif
	
	else 
	
		max_limiter = max(q_cellface_y(iX,iY-1,1),q_cellface_y(iX,iY,1),q_half_step) + source_term
		min_limiter = min(q_cellface_y(iX,iY-1,1),q_cellface_y(iX,iY,1),q_half_step) + source_term
		
		if (q_cellface_y(iX,iY,2) > max_limiter) then
			
			q_cellface_y(iX,iY,2) = max_limiter
		
		endif
		
		if (q_cellface_y(iX,iY,2) < min_limiter) then
			
			q_cellface_y(iX,iY,2) = min_limiter
		
		endif
		
		
	endif


end subroutine flux_limitor_y2

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



	
end program TestCABARET2D





