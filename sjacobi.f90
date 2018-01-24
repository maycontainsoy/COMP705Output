! Solve Laplace equation using Jacobi iteration method
! Kadin Tseng, Boston University, November 1999

PROGRAM Jacobi
USE jacobi_module
REAL(real8), DIMENSION(:,:), POINTER :: c, n, e, w, s
INTEGER :: i,iRow
CHARACTER(len=10) :: arg, charIter

K = 0                        ! thread #
P = 1                        ! # of threads
!write(*,"('Enter size of interior points, m : ')",advance='NO')
!read(*,*)m
DO i = 1,iargc()
  CALL GETARG(i,arg)
  IF ( i .EQ. 1 ) THEN
    READ (arg,*) m
  ELSE IF ( i .EQ. 2 ) THEN
    READ (arg,*) iterSave
  ELSE
    testName = arg 
    PRINT *, testName
  END IF
END DO

mp = m/P

CALL cpu_time(start_time)    ! start timer, measured in seconds

ALLOCATE ( vnew(m,m), v(0:m+1,0:m+1) )  ! mem for vnew, v

c => v(1:m  ,1:m  )          ! i  ,j    Current/Central for 1<=i<=m;1<=j<=m
n => v(1:m  ,2:m+1)          ! i  ,j+1  North (of Current)
e => v(2:m+1,1:m  )          ! i+1,j    East  (of Current)
w => v(0:m-1,1:m  )          ! i-1,j    West  (of Current)
s => v(1:m  ,0:m-1)          ! i  ,j-1  South (of Current)

CALL bc(v, m, mp, K, P)      ! set up boundary values

DO WHILE (gdel > TOL)        ! iterate until error below threshold
  iter = iter + 1            ! increment iteration counter
  IF(iter > MAXSTEPS) THEN
    WRITE(*,*)'Iteration terminated (exceeds', MAXSTEP,')'
    STOP                     ! nonconvergent solution
  ENDIF
  vnew = ( n + e + w + s )*0.25 ! new solution, Eq. 3
  gdel = SUM(DABS(vnew-c))   ! find local max error
  IF(MOD(iter,INCREMENT)==0) WRITE(*,"('iter,gdel:',i6,e12.4)")iter,gdel

  IF(MOD(iter,iterSave)==0) THEN
    PRINT *, 'Mod check for iterSave'
    WRITE(charIter,'(I5.5)') iter
    PRINT *, 'Character Iteration: ', charIter
    outfile = 'plots' // '.' // TRIM(testName) // '.' // TRIM(charIter)
    OPEN(unit=39,file=outfile,form='formatted',status='unknown')
    DO iRow=0,m+1
!      WRITE(39,*) v(:,:)
      WRITE(39,"(e13.4)") v(:,iRow)
    END DO
!    WRITE(39,"(6e13.4)") v
    CLOSE(39)
  END IF
  c = vnew                   ! update interior v
ENDDO

CALL CPU_TIME(end_time)      ! stop timer
PRINT *,'Total cpu time =',end_time - start_time,' x 1'
PRINT *,'Stopped at iteration =',iter
PRINT *,'The maximum error =',gdel

!OPEN(unit=40, file=outfile, form='formatted', status='unknown')
!WRITE(40,"(3i5)")m+2,m+2,P
!close(40)

PRINT *, 'Writes to original plots file'

!WRITE(outfile, "('output.',a)") testName
outfile = TRIM(outfile) // '.' // TRIM(testName)
!WRITE(outfile, '(a,a)') "plots",testName
PRINT *, 'Plot output file:', outfile
OPEN(unit=41, file=outfile, form='formatted', status='unknown')
WRITE(41,"(6e13.4)")v
close(41)

PRINT *, 'Writes to plots.exe file'

DEALLOCATE (vnew, v)

END PROGRAM Jacobi
