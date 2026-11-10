!Copyright 2008-2011 SEISCOPE project, All rights reserved.
!Copyright 2013-2021 SEISCOPEII project, All rights reserved.
subroutine sub_read_freqdom_freq_list(pbdir)

  implicit none


#include "pbdirect.h"
  include "mpif.h"
#include "common.h"


  TYPE (pbdirect) :: pbdir
  INTEGER :: i

  OPEN(12,file='freq_management')
  READ(12,*)pbdir%nw
  ALLOCATE(pbdir%freq(pbdir%nw))
  ALLOCATE(pbdir%freqw(pbdir%nw))
  READ(12,*)(pbdir%freq(i),i=1,pbdir%nw)
  CLOSE(12)
  DO i=1,pbdir%nw
     pbdir%freqw(i)=pbdir%freq(i)*2.*pi
  END DO

end subroutine sub_read_freqdom_freq_list
