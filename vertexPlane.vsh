#version 330

uniform mat4 vpMatrix;
layout(location=0) in vec3 vertex_modelspace;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace,1);
}
