SUBROUTINE make_iso_surf(volume,thresh,x,y,z,vertices,nv)

USE iso_surf_lut
IMPLICIT NONE

INTEGER,INTENT(IN)::x,y,z
REAL,INTENT(IN)::thresh
REAL,INTENT(IN),DIMENSION(x,y,z)::volume

INTEGER,INTENT(OUT)::nv
REAL,INTENT(OUT),DIMENSION(3,3,x*y*z*4)::vertices

INTEGER::i,i2,j,j2,k,k2,l,vsum
REAL,DIMENSION(3)::cellshape
REAL,DIMENSION(3,0:11)::tempvert

cellshape=(/x,y,z/)
nv=1
DO i=1,x
    !care for PBC
    IF(i==x)THEN
        i2=1
    ELSE
        i2=i+1
    END IF
    DO j=1,y
        IF(j==y)THEN
            j2=1
        ELSE
            j2=j+1
        END IF
        DO k=1,z
            IF(k==z)THEN
                k2=1
            ELSE
                k2=k+1
            END IF

            !determine the cube-vertices inside of volume
            vsum=0
            IF(volume(i ,j ,k )<thresh) vsum=vsum+1
            IF(volume(i ,j ,k2)<thresh) vsum=vsum+2
            IF(volume(i ,j2,k )<thresh) vsum=vsum+4
            IF(volume(i ,j2,k2)<thresh) vsum=vsum+8
            IF(volume(i2,j ,k )<thresh) vsum=vsum+16
            IF(volume(i2,j ,k2)<thresh) vsum=vsum+32
            IF(volume(i2,j2,k )<thresh) vsum=vsum+64
            IF(volume(i2,j2,k2)<thresh) vsum=vsum+128

            IF (0<vsum .and. vsum<255) THEN

                !determine the interpolated edge-vertices for intersection
                tempvert=0
                IF(IAND(edge_lut(vsum),   1)==1)THEN
                    tempvert(:, 0)=vert_off(:,0)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,0)-vert_off(:,1))
                END IF
                IF(IAND(edge_lut(vsum),   2)==1)THEN
                    tempvert(:, 1)=vert_off(:,1)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,1)-vert_off(:,5))
                END IF
                IF(IAND(edge_lut(vsum),   4)==1)THEN
                    tempvert(:, 2)=vert_off(:,5)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,5)-vert_off(:,4))
                END IF
                IF(IAND(edge_lut(vsum),   8)==1)THEN
                    tempvert(:, 3)=vert_off(:,0)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,0)-vert_off(:,4))
                END IF
                IF(IAND(edge_lut(vsum),  16)==1)THEN
                    tempvert(:, 4)=vert_off(:,1)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,1)-vert_off(:,3))
                END IF
                IF(IAND(edge_lut(vsum),  32)==1)THEN
                    tempvert(:, 5)=vert_off(:,5)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,5)-vert_off(:,7))
                END IF
                IF(IAND(edge_lut(vsum),  64)==1)THEN
                    tempvert(:, 6)=vert_off(:,4)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,4)-vert_off(:,6))
                END IF
                IF(IAND(edge_lut(vsum), 128)==1)THEN
                    tempvert(:, 7)=vert_off(:,0)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,0)-vert_off(:,2))
                END IF
                IF(IAND(edge_lut(vsum), 256)==1)THEN
                    tempvert(:, 8)=vert_off(:,2)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,2)-vert_off(:,3))
                END IF
                IF(IAND(edge_lut(vsum), 512)==1)THEN
                    tempvert(:, 9)=vert_off(:,3)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,3)-vert_off(:,7))
                END IF
                IF(IAND(edge_lut(vsum),1024)==1)THEN
                    tempvert(:,10)=vert_off(:,6)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,6)-vert_off(:,7))
                END IF
                IF(IAND(edge_lut(vsum),2048)==1)THEN
                    tempvert(:,11)=vert_off(:,2)+&
                        ((thresh-volume(i,j,k))/(volume(i,j,k2)-volume(i,j,k)))*(vert_off(:,2)-vert_off(:,6))
                END IF

                !create triangles
                DO l=1,nv_lut(vsum)
                    vertices(:,1,nv)=((/i-1,j-1,k-1/)+tempvert(:,tri_lut(1,l,vsum)))/(/x,y,z/)
                    vertices(:,2,nv)=((/i-1,j-1,k-1/)+tempvert(:,tri_lut(2,l,vsum)))/(/x,y,z/)
                    vertices(:,3,nv)=((/i-1,j-1,k-1/)+tempvert(:,tri_lut(3,l,vsum)))/(/x,y,z/)
                    nv=nv+1
                END DO
            END IF

        END DO
    END DO
END DO
nv=nv-1

END SUBROUTINE make_iso_surf
