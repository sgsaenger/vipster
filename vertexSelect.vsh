#version 130

//! [0]
uniform mat4 mvpMatrix;
in vec3 vertex_modelspace;

void main(void)
{
    //standard vertex positioning:
    gl_Position = mvpMatrix * vec4(vertex_modelspace,1);
}
//! [0]
