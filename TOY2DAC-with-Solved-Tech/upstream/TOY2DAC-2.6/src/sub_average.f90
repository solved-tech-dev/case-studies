!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
!
! subroutine subbaverage
! interpolate buoyancy at intermediate grid points
! arithmetic averaging is used here.
! b(0:n1e+1,0:n2e+1): buoyancy array on the reference grid
! bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m: interpolated values
! bpm=b(i1+1/2,i2-1/2), bpp=b(i1+1/2,i2+1/2), ...

subroutine subbaverage(b,n1e,n2e,i1,i2,bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m)
integer n1e,n2e,i1,i2
real b(0:n1e+1,0:n2e+1)
real bpm,bmp,bpp,bmm,b00,bp0,bm0,b0p,b0m

bpm=0.25*(b(i1,i2)+b(i1+1,i2)+b(i1,i2-1)+b(i1+1,i2-1))
bmp=0.25*(b(i1,i2)+b(i1-1,i2)+b(i1,i2+1)+b(i1-1,i2+1))
bpp=0.25*(b(i1,i2)+b(i1+1,i2)+b(i1,i2+1)+b(i1+1,i2+1))
bmm=0.25*(b(i1,i2)+b(i1-1,i2)+b(i1,i2-1)+b(i1-1,i2-1))
b00=b(i1,i2)
bp0=0.5*(b(i1,i2)+b(i1+1,i2))
bm0=0.5*(b(i1,i2)+b(i1-1,i2))
b0p=0.5*(b(i1,i2)+b(i1,i2+1))
b0m=0.5*(b(i1,i2)+b(i1,i2-1))

return
end
