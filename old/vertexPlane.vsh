#version 330

uniform mat4 vpMatrix;
layout(location=0) in vec3 vertex_modelspace;
layout(location=1) in vec2 vertex_UV;
layout(location=2) in vec3 offset;

out vec2 UV;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex_modelspace+offset,1);
    UV=vertex_UV;
}
