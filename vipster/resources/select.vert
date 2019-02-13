layout(location = 0) in vec3 vertex;
layout(location = 1) in vec3 position;
layout(location = 2) in float vert_scale;
layout(location = 4) in uint hide;

layout(std140, row_major) uniform viewMat{
    mat4 vpMatrix;
    mat4 rMatrix;
};
uniform mat3 pos_scale;
uniform float scale_fac;
uniform vec3 offset;

out vec4 color_input;
flat out uint frag_discard;

void main(void)
{
    frag_discard = hide;
    if(hide==uint(0)){
        gl_Position = vpMatrix * vec4(vertex * vert_scale * scale_fac +
                                      position * pos_scale +
                                      offset, 1);
        float red=float(gl_InstanceID&0xFF)/255.;
        float green=float((gl_InstanceID&0xFF00)>>8)/255.f;
        float blue=float((gl_InstanceID&0xFF0000)>>16)/255.f;
        color_input = vec4(red,green,blue,1);
    }
}
