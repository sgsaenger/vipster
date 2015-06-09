SUBROUTINE make_vol_gradient(volume,x,y,z,gradient)

IMPLICIT NONE
INTEGER,INTENT(IN)::x,y,z
REAL,INTENT(IN),DIMENSION(x,y,z)::volume

REAL,INTENT(OUT),DIMENSION(3,x,y,z)::gradient

INTEGER::i,j,k
INTEGER::il,ih,jl,jh,kl,kh
!REAL::gx,gy,gz,norm

DO i=1,x
    !care for PBC
    IF(i==1)THEN
        il=x
        ih=2
    ELSE IF(i==x)THEN
        il = x-1
        ih = 1
    ELSE
        il = i-1
        ih = i+1
    END IF
    DO j=1,y
        IF(j==1)THEN
            jl=y
            jh=2
        ELSE IF(j==y)THEN
            jl = y-1
            jh = 1
        ELSE
            jl = j-1
            jh = j+1
        END IF
        DO k=1,z
            IF(k==1)THEN
                kl=z
                kh=2
            ELSE IF(k==z)THEN
                kl = z-1
                kh = 1
            ELSE
                kl = k-1
                kh = k+1
            END IF

            !calc normalized gradient
            gradient(1,i,j,k)=volume(il,j,k)-volume(ih,j,k)
            gradient(2,i,j,k)=volume(i,jl,k)-volume(i,jh,k)
            gradient(3,i,j,k)=volume(i,j,kl)-volume(i,j,kh)
        END DO
    END DO
END DO

END SUBROUTINE
