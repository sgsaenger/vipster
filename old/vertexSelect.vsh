in vec3 vertex_modelspace;
#if __VERSION__ < 330
uniform vec3 position_modelspace;
uniform float scale_modelspace;
uniform vec4 in_color;
#else
in vec3 position_modelspace;
in float scale_modelspace;
#endif

uniform mat4 vpMatrix;

out vec4 MaterialDiffuseColor;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace*scale_modelspace+position_modelspace,1);

#if __VERSION__ < 330
    MaterialDiffuseColor = in_color;
#else
    float red=float(gl_InstanceID&0xFF)/255.;
    float green=float((gl_InstanceID&0xFF00)>>8)/255.;
    float blue=float((gl_InstanceID&0xFF0000)>>16)/255.;
    MaterialDiffuseColor = vec4(red,green,blue,1);
#endif
}
