#version 130

//! [0]
uniform vec4 color;

out vec3 cellColor;

void main(void)
{
    cellColor = color.xyz;
}
//! [0]
