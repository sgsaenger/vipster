layout(location = 0) in vec3 vertex;

uniform mat3 pos_scale;
uniform vec3 offset;
uniform vec4 color;
layout(std140, row_major) uniform viewMat{
    mat4 vpMatrix;
    mat4 rMatrix;
};

out vec4 color_input;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex * pos_scale + offset, 1);
//    gl_Position = vpMatrix * vec4(vertex * pos_scale + offset, 1);
    color_input = color;
}
