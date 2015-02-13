#version 330

uniform mat4 mvpMatrix;
uniform mat4 mvMatrix;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;

layout(location=2) in vec3 vertex_modelspace;
//normals == vertices

void main(void)
{
    //standard vertex positioning:
    gl_Position = mvpMatrix * vec4(vertex_modelspace,1);

    //transformations for calculation of lighting:
    vec3 vertex_cameraspace=(mvMatrix*vec4(vertex_modelspace,1)).xyz;

    //from vertex to camera (in cameraspace always at origin)
    EyeDirection_cameraspace = vec3(0,0,0) - vertex_cameraspace;

    //from vertex to light source
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;

    //Normals in camera space
    normals_cameraspace = (mvMatrix * vec4(vertex_modelspace,0)).xyz;
}
