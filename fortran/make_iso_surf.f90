SUBROUTINE make_iso_surf(volume,thresh,x,y,z,gradient,vertices,nv)

USE iso_surf_lut
IMPLICIT NONE

INTEGER,INTENT(IN)::x,y,z
REAL,INTENT(IN)::thresh
REAL,INTENT(IN),DIMENSION(x,y,z)::volume
REAL,INTENT(IN),DIMENSION(3,x,y,z)::gradient

INTEGER,INTENT(OUT)::nv
REAL,INTENT(OUT),DIMENSION(3,6,x*y*z*4)::vertices

INTEGER::i,i2,j,j2,k,k2,l,vsum
REAL::ratio
REAL,DIMENSION(3)::cellshape
REAL,DIMENSION(3,0:11)::tempvert,tempnorm

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
                tempnorm=0
                IF(IAND(edge_lut(vsum),   1)/=0)THEN
                    ratio=((thresh-volume(i ,j ,k ))/(volume(i ,j ,k2)-volume(i ,j ,k )))
                    tempvert(:, 0)=vert_off(:,0)+ratio*(vert_off(:,1)-vert_off(:,0))
                    tempnorm(:, 0)=gradient(:,i ,j ,k )+ratio*(gradient(:,i ,j ,k2)-gradient(:,i ,j ,k ))
                END IF
                IF(IAND(edge_lut(vsum),   2)/=0)THEN
                    ratio=((thresh-volume(i ,j ,k2))/(volume(i2,j ,k2)-volume(i ,j ,k2)))
                    tempvert(:, 1)=vert_off(:,1)+ratio*(vert_off(:,5)-vert_off(:,1))
                    tempnorm(:, 1)=gradient(:,i ,j ,k2)+ratio*(gradient(:,i2,j ,k2)-gradient(:,i ,j ,k2))
                END IF
                IF(IAND(edge_lut(vsum),   4)/=0)THEN
                    ratio=((thresh-volume(i2,j ,k2))/(volume(i2,j ,k )-volume(i2,j ,k2)))
                    tempvert(:, 2)=vert_off(:,5)+ratio*(vert_off(:,4)-vert_off(:,5))
                    tempnorm(:, 2)=gradient(:,i2,j ,k2)+ratio*(gradient(:,i2,j ,k )-gradient(:,i2,j ,k2))
                END IF
                IF(IAND(edge_lut(vsum),   8)/=0)THEN
                    ratio=((thresh-volume(i ,j ,k ))/(volume(i2,j ,k )-volume(i ,j ,k )))
                    tempvert(:, 3)=vert_off(:,0)+ratio*(vert_off(:,4)-vert_off(:,0))
                    tempnorm(:, 3)=gradient(:,i ,j ,k )+ratio*(gradient(:,i2,j ,k )-gradient(:,i ,j ,k ))
                END IF
                IF(IAND(edge_lut(vsum),  16)/=0)THEN
                    ratio=((thresh-volume(i ,j ,k2))/(volume(i ,j2,k2)-volume(i ,j ,k2)))
                    tempvert(:, 4)=vert_off(:,1)+ratio*(vert_off(:,3)-vert_off(:,1))
                    tempnorm(:, 4)=gradient(:,i ,j ,k2)+ratio*(gradient(:,i ,j2,k2)-gradient(:,i ,j ,k2))
                END IF
                IF(IAND(edge_lut(vsum),  32)/=0)THEN
                    ratio=((thresh-volume(i2,j ,k2))/(volume(i2,j2,k2)-volume(i2,j ,k2)))
                    tempvert(:, 5)=vert_off(:,5)+ratio*(vert_off(:,7)-vert_off(:,5))
                    tempnorm(:, 5)=gradient(:,i2,j ,k2)+ratio*(gradient(:,i2,j2,k2)-gradient(:,i2,j ,k2))
                END IF
                IF(IAND(edge_lut(vsum),  64)/=0)THEN
                    ratio=((thresh-volume(i2,j ,k ))/(volume(i2,j2,k )-volume(i2,j ,k )))
                    tempvert(:, 6)=vert_off(:,4)+ratio*(vert_off(:,6)-vert_off(:,4))
                    tempnorm(:, 6)=gradient(:,i2,j ,k )+ratio*(gradient(:,i2,j2,k )-gradient(:,i2,j ,k ))
                END IF
                IF(IAND(edge_lut(vsum), 128)/=0)THEN
                    ratio=((thresh-volume(i ,j ,k ))/(volume(i ,j2,k )-volume(i ,j ,k )))
                    tempvert(:, 7)=vert_off(:,0)+ratio*(vert_off(:,2)-vert_off(:,0))
                    tempnorm(:, 7)=gradient(:,i ,j ,k )+ratio*(gradient(:,i ,j2,k )-gradient(:,i ,j ,k ))
                END IF
                IF(IAND(edge_lut(vsum), 256)/=0)THEN
                    ratio=((thresh-volume(i ,j2,k ))/(volume(i ,j2,k2)-volume(i ,j2,k )))
                    tempvert(:, 8)=vert_off(:,2)+ratio*(vert_off(:,3)-vert_off(:,2))
                    tempnorm(:, 8)=gradient(:,i ,j2,k )+ratio*(gradient(:,i ,j2,k2)-gradient(:,i ,j2,k ))
                END IF
                IF(IAND(edge_lut(vsum), 512)/=0)THEN
                    ratio=((thresh-volume(i ,j2,k2))/(volume(i2,j2,k2)-volume(i ,j2,k2)))
                    tempvert(:, 9)=vert_off(:,3)+ratio*(vert_off(:,7)-vert_off(:,3))
                    tempnorm(:, 9)=gradient(:,i ,j2,k2)+ratio*(gradient(:,i2,j2,k2)-gradient(:,i ,j2,k2))
                END IF
                IF(IAND(edge_lut(vsum),1024)/=0)THEN
                    ratio=((thresh-volume(i2,j2,k ))/(volume(i2,j2,k2)-volume(i2,j2,k )))
                    tempvert(:,10)=vert_off(:,6)+ratio*(vert_off(:,7)-vert_off(:,6))
                    tempnorm(:,10)=gradient(:,i2,j2,k )+ratio*(gradient(:,i2,j2,k2)-gradient(:,i2,j2,k ))
                END IF
                IF(IAND(edge_lut(vsum),2048)/=0)THEN
                    ratio=((thresh-volume(i ,j2,k ))/(volume(i2,j2,k )-volume(i ,j2,k )))
                    tempvert(:,11)=vert_off(:,2)+ratio*(vert_off(:,6)-vert_off(:,2))
                    tempnorm(:,11)=gradient(:,i ,j2,k )+ratio*(gradient(:,i2,j2,k )-gradient(:,i ,j2,k ))
                END IF

                !create triangles
                DO l=1,nv_lut(vsum)
                    vertices(:,1,nv)=((/i-1,j-1,k-1/)+tempvert(:,tri_lut(1,l,vsum)))/(/x,y,z/)
                    vertices(:,2,nv)=tempnorm(:,tri_lut(1,l,vsum))
                    vertices(:,3,nv)=((/i-1,j-1,k-1/)+tempvert(:,tri_lut(2,l,vsum)))/(/x,y,z/)
                    vertices(:,4,nv)=tempnorm(:,tri_lut(2,l,vsum))
                    vertices(:,5,nv)=((/i-1,j-1,k-1/)+tempvert(:,tri_lut(3,l,vsum)))/(/x,y,z/)
                    vertices(:,6,nv)=tempnorm(:,tri_lut(3,l,vsum))
                    nv=nv+1
                END DO
            END IF

        END DO
    END DO
END DO
nv=nv-1

END SUBROUTINE make_iso_surf
