#version 330 core

layout(location = 0) in vec3 vertex_modelspace;
uniform mat4 vpMatrix;
uniform vec3 offset;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex_modelspace+offset,1);
}

