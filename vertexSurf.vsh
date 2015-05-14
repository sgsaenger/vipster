#version 330

uniform mat4 vpMatrix;
//uniform mat3 cellVec;
uniform vec3 shape;
uniform isamplerBuffer vol_val;

layout(location=0) in vec3 offset;

out int vertices;
//out vec3 coord;
out int vID;

void main(void)
{
    vID = gl_VertexID;
    //coord.z = (0.5+mod(gl_VertexID,shape.z))/shape.z;
    //coord.y = (0.5+mod(gl_VertexID/int(shape.z),shape.y))/shape.y;
    //coord.x = (0.5+mod(gl_VertexID/int(shape.z*shape.y),shape.x))/shape.x;
    //gl_Position = vpMatrix * vec4(cellVec*coord+offset,1);
    gl_Position = vec4(offset,0);

    vertices=texelFetch(vol_val,gl_VertexID).r;
}
