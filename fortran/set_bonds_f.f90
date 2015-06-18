SUBROUTINE set_bonds_f(nat,coord,cut,cut_fac,off,k,at1,at2,dist)
IMPLICIT NONE
INTEGER :: nat,i,j,k
!f2py intent(in) nat
!f2py intent(out) k
REAL,INTENT(IN)::cut_fac
REAL,INTENT(IN),DIMENSION(2,3)::off
REAL,INTENT(IN),DIMENSION(0:nat-1,3) :: coord
REAL,INTENT(IN),DIMENSION(0:nat-1) :: cut
INTEGER,INTENT(OUT),DIMENSION(nat*nat) :: at1,at2
REAL,INTENT(OUT),DIMENSION(nat*nat) :: dist
REAL,DIMENSION(3) :: dist_v
REAL :: effcut,dist_n
!add one in order to use it as an index in fortran
k=1
DO i=0,nat-1
        DO j=0,nat-1
                IF ((j<i).AND.ALL(off==0.)) CYCLE
                IF (i==j) CYCLE
                dist_v = coord(i,:)+off(1,:)-coord(j,:)-off(2,:)
                effcut=(cut(i)+cut(j))*cut_fac
                IF(ANY(dist_v>effcut))CYCLE
                dist_n=DOT_PRODUCT(dist_v,dist_v)
                IF((0.57<dist_n).AND.(dist_n<effcut**2))THEN
                        at1(k)=i
                        at2(k)=j
                        dist(k)=sqrt(dist_n)
                        k=k+1
                END IF
        END DO
END DO
!remove one in order to use it as an index in python
k=k-1

END SUBROUTINE set_bonds_f
