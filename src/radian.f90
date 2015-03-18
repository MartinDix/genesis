!------ radian ------
!  convert degrees to radians

real function radian(ang)
	
    implicit none
	
    real pi
    parameter(pi=3.1415926)
    real ang
        
    radian = ang*pi/180.

    return
end function radian
        
!-------------------------------------------        
