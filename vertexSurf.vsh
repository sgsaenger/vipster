#version 330

uniform mat4 vpMatrix;
uniform mat4 rMatrix;
uniform mat3 cellVec;
uniform vec3 volOff;

layout(location=0) in vec3 vertex_cellspace;
layout(location=1) in vec3 normal_cellspace;
layout(location=2) in vec3 offset;

out vec4 MaterialDiffuseColor;
out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * vec4(vertex_cellspace*cellVec+offset+volOff,1);
    //fixed color for now
    MaterialDiffuseColor = vec4(0.8,0.1,0.1,1);

    //prepare normals:
    vec3 vertex_cameraspace=(rMatrix*vec4(vertex_cellspace*cellVec,1)).xyz;
    EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;
    normals_cameraspace = (rMatrix * vec4(normal_cellspace,0)).xyz;
}
