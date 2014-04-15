!------------------------------------------------------------
!
	real function modulus(dd)
	
	implicit none
	
	real,intent(in)		:: dd
	
	modulus = sqrt(dd**2)
	
	return
	
	end function modulus
	
	
