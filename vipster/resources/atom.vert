//layout(location = 0) in vec3 vertex;
layout(location = 1) in vec3 position;
layout(location = 2) in float vert_scale;
layout(location = 3) in vec4 color;
layout(location = 4) in uint hide;

layout(std140, row_major) uniform viewMat{
    mat4 rvpMatrix;
    mat4 pMatrix;
    mat4 vMatrix;
    mat4 rMatrix;
};
uniform mat3 pos_scale;
uniform float scale_fac;
uniform vec3 offset;
uniform vec4 viewport;

out vec4 MaterialDiffuseColor;
flat out uint frag_discard;
out float radius;
out vec2 center;

void main(void)
{
    frag_discard = hide;
    if(hide==uint(0)){
        gl_Position = rvpMatrix * vec4(position * pos_scale + offset, 1);
        center = gl_Position.xy;
        // take scaling from vMatrix
        // (assumed to contain only scaling and translation)
        radius = vMatrix[0][0] * vert_scale * scale_fac;
        gl_PointSize = radius * min(viewport.z, viewport.w);
        MaterialDiffuseColor = color;
    }
}

