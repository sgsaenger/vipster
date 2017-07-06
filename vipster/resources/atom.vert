layout(location = 0) in vec3 vertex_modelspace;
layout(location = 1) in vec3 position_modelspace;
layout(location = 2) in float scale_modelspace;
layout(location = 3) in vec4 color_input;

uniform mat4 vpMatrix;
uniform mat4 rMatrix;
uniform float atom_fac;
uniform vec3 offset;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out vec4 MaterialDiffuseColor;

void main(void)
{
    gl_Position = vpMatrix *vec4(vertex_modelspace*scale_modelspace*atom_fac+position_modelspace+offset,1);
    vec3 vertex_cameraspace=(rMatrix*vec4(vertex_modelspace,1)).xyz;
    EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;
    normals_cameraspace = (rMatrix*vec4(vertex_modelspace,0)).xyz;
    MaterialDiffuseColor = color_input;
}

