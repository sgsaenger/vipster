SUBROUTINE make_iso_surf(volume,thresh,x,y,z,vertices,nv)
USE iso_surf_lut
IMPLICIT NONE
INTEGER,INTENT(IN)::x,y,z
REAL,INTENT(IN)::thresh
REAL,INTENT(IN),DIMENSION(x,y,z)::volume

INTEGER,INTENT(OUT)::nv
REAL,INTENT(OUT),DIMENSION(3,3,x*y*z*4)::vertices

INTEGER::i,j,k,l
INTEGER,DIMENSION(x,y,z)::v1,v2,v4,v8,v16,v32,v64,v128,vsum
IF(thresh>0)THEN
    v1=volume>thresh    !Bit 0: 1
ELSE IF(thresh<0)THEN
    v1=volume<thresh    !Bit 0: 1
END IF
v2 = cshift(v1,1,3)*2  !Bit 1: 2 z
v4 = cshift(v1,1,2)*4  !Bit 2: 4 y
v8 = cshift(v2,1,2)*4  !Bit 3: 8 = 2*4 yz
v16 = cshift(v1,1,1)*16    !Bit 4: 16 x
v32 = cshift(v2,1,1)*16    !Bit 5: 32 = 2*16 xz
v64 = cshift(v4,1,1)*16    !Bit 6: 64 = 4*16 xy
v128 = cshift(v8,1,1)*16   !Bit 7: 128 = 8*16 xyz

vsum = v1+v2+v4+v8+v16+v32+v64+v128

nv = 1
vertices=0.

DO i=1,x
    DO j=1,y
        DO k=1,z
            IF (0<vsum(i,j,k) .and. vsum(i,j,k)<255) THEN
                DO l=1,nv_lut(vsum(i,j,k))
                    vertices(:,1,nv)=((/i-1,j-1,k-1/)+edge_off(:,edge_lut(1,l,vsum(i,j,k))))/(/x,y,z/)
                    vertices(:,2,nv)=((/i-1,j-1,k-1/)+edge_off(:,edge_lut(2,l,vsum(i,j,k))))/(/x,y,z/)
                    vertices(:,3,nv)=((/i-1,j-1,k-1/)+edge_off(:,edge_lut(3,l,vsum(i,j,k))))/(/x,y,z/)
                    nv=nv+1
                END DO
            END IF
        END DO
    END DO
END DO
nv = nv-1

END SUBROUTINE make_iso_surf
