layout(location = 0) in vec3 vertex;

uniform vec3 offset;
uniform vec4 color;

out vec4 color_input;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex+offset, 1);
    color_input = color;
}
