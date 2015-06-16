uniform mat4 vpMatrix;
in vec3 vertex_modelspace;
in vec2 vertex_UV;
#if __VERSION__ < 330
uniform vec3 offset;
#else
in vec3 offset;
#endif

out vec2 UV;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex_modelspace+offset,1);
    UV=vertex_UV;
}
