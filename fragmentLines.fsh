#version 330

//! [0]
uniform vec4 color;

out vec4 cellColor;

void main(void)
{
    cellColor = color;
}
//! [0]
