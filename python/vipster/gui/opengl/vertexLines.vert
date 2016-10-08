uniform mat4 vpMatrix;
in vec3 vertex_modelspace;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace,1);
}
