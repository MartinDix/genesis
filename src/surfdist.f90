!------ surfdist -------
!great circle distance between lat/lon pairs
!

real function surfdist(lat1,lon1,lat2,lon2)

  implicit none

  real radius
  parameter (radius = 6371220)
  real lon1,lat1,lon2,lat2
  real cosang
  real ang
  real rad
  real radian

  cosang = cos(radian(lon2)-radian(lon1))
  cosang = cosang * cos(radian(lat1)) * cos(radian(lat2))
  cosang = cosang + sin(radian(lat1))*sin(radian(lat2))

  ang = acos(cosang)
  if (ang .lt. 0 ) stop',Neg center angle!!'
  surfdist = ang*radius

  return
end function surfdist
