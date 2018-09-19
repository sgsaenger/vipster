layout(location = 0) in vec3 vertex;
layout(location = 1) in vec3 normal;
layout(location = 2) in vec2 vert_UV;

uniform mat3 pos_scale;
uniform vec3 offset;
layout(std140, row_major) uniform viewMat{
    mat4 vpMatrix;
    mat4 rMatrix;
};

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out vec2 UV;

void main(void)
{
    gl_Position = vpMatrix * vec4(vertex * pos_scale + offset, 1);
    vec3 vertex_cameraspace = (rMatrix * vec4(vertex, 1)).xyz;
    EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;
    normals_cameraspace = (rMatrix * vec4(normal, 0)).xyz;
    UV = vert_UV;
}
