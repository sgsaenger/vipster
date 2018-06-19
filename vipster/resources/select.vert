layout(location = 0) in vec3 vertex_modelspace;
layout(location = 1) in vec3 position_modelspace;
layout(location = 2) in float scale_modelspace;

layout(std140, row_major) uniform viewMat{
    mat4 vpMatrix;
    mat4 rMatrix;
};
uniform mat3 position_scale;
uniform float atom_fac;
uniform vec3 offset;

out vec4 color_input;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex_modelspace * scale_modelspace * atom_fac +
                                  position_modelspace * position_scale +
                                  offset, 1);
    float red=float(gl_InstanceID&0xFF)/255.;
    float green=float((gl_InstanceID&0xFF00)>>8)/255.f;
    float blue=float((gl_InstanceID&0xFF0000)>>16)/255.f;
    color_input = vec4(red,green,blue,1);
}
