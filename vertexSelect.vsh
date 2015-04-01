#version 330

layout(location=0) in vec3 vertex_modelspace;
layout(location=1) in vec3 position_modelspace;
layout(location=2) in float scale_modelspace;

uniform mat4 vpMatrix;

out vec4 MaterialDiffuseColor;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace*scale_modelspace+position_modelspace,1);

    //MaterialDiffuseColor = vec4(float((gl_InstanceID&16711680)>>16)/255.,float((gl_InstanceID&65280)>>8)&255.,float(gl_InstanceID&255)/255.,1);
    float red=float(gl_InstanceID&0xFF)/255.;
    float green=float((gl_InstanceID&0xFF00)>>8)/255.;
    float blue=float((gl_InstanceID&0xFF0000)>>16)/255.;
    MaterialDiffuseColor = vec4(red,green,blue,1);
}
