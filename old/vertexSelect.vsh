#version 330

layout(location=0) in vec3 vertex_modelspace;
uniform vec3 position_modelspace;
uniform float scale_modelspace;
uniform vec4 in_color;

uniform mat4 vpMatrix;

out vec4 MaterialDiffuseColor;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace*scale_modelspace+position_modelspace,1);

    MaterialDiffuseColor = in_color;
}
