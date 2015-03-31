SUBROUTINE set_bonds_f(nat,ci,ni,off,k,co1,co2,no1,no2)
IMPLICIT NONE
INTEGER :: nat,i,j,k
!f2py intent(in) nat
!f2py intent(out) k
REAL,INTENT(IN),DIMENSION(2,3) :: off
REAL,INTENT(IN),DIMENSION(nat,3) :: ci
CHARACTER,INTENT(IN),DIMENSION(nat) :: ni
REAL,INTENT(OUT),DIMENSION((nat*(nat-1))/2,3) :: co1,co2
CHARACTER,INTENT(OUT),DIMENSION((nat*(nat-1))/2) :: no1,no2
REAL,DIMENSION(3) :: dist_v
REAL :: dist_n
!add one in order to use it as an index in fortran
k=1
DO i=1,nat
        DO j=i+1,nat
                dist_v = ci(i,:)+off(1,:)-ci(j,:)-off(2,:)
                IF((dist_v(1)>3.5).OR.(dist_v(2)>3.5).OR.(dist_v(3)>3.5)) CYCLE
                dist_n = dot_product(dist_v,dist_v)
                IF((ni(i)/='H').AND.(ni(j)/='H'))THEN
                        IF((0.57<dist_n).AND.(dist_n<12.25))THEN
                                co1(k,:)=ci(i,:)+off(1,:)
                                co2(k,:)=ci(j,:)+off(2,:)
                                no1(k)=ni(i)
                                no2(k)=ni(j)
                                k=k+1
                        END IF
                ELSE
                        IF((0.57<dist_n).AND.(dist_n<5.15))THEN
                                co1(k,:)=ci(i,:)+off(1,:)
                                co2(k,:)=ci(j,:)+off(2,:)
                                no1(k)=ni(i)
                                no2(k)=ni(j)
                                k=k+1
                        END IF
                END IF
        END DO
END DO
!remove one in order to use it as an index in python
k=k-1
END SUBROUTINE set_bonds_f
