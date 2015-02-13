#version 330

layout(location=0) in vec3 position_modelspace;
layout(location=1) in float scale_modelspace;
layout(location=2) in vec3 vertex_modelspace;
layout(location=3) in vec4 color_input;
//normals == vertices

uniform mat4 vpMatrix;
uniform mat4 vMatrix;
uniform mat4 rMatrix;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out vec4 MaterialDiffuseColor;


void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_modelspace*scale_modelspace+position_modelspace,1);

    //transformations for calculation of lighting:
    vec3 vertex_cameraspace=(rMatrix*vec4(vertex_modelspace,1)).xyz;

    //from vertex to camera (in cameraspace always at origin)
    EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;

    //from vertex to light source
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;

    //Normals in camera space
    normals_cameraspace = (rMatrix * vec4(vertex_modelspace,0)).xyz;

    // Pass instanced color to fragmentshader
    MaterialDiffuseColor = color_input;
}
