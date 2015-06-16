#version 330

layout(location=0) in vec3 vertex_modelspace;
uniform vec3 position_modelspace;
uniform float scale_modelspace;
uniform int InstanceID;

uniform mat4 vpMatrix;

out vec4 MaterialDiffuseColor;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace*scale_modelspace+position_modelspace,1);

    float red=float(InstanceID&0xFF)/255.;
    float green=float((InstanceID&0xFF00)>>8)/255.;
    float blue=float((InstanceID&0xFF0000)>>16)/255.;
    MaterialDiffuseColor = vec4(red,green,blue,1);
}
