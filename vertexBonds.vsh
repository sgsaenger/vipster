#version 130

//! [0]
uniform mat4 mvpMatrix;
uniform mat4 vMatrix;
uniform mat4 mMatrix;
uniform vec3 LightPosition_cameraspace;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out float vertex_side;

in vec3 vertex_modelspace;
in vec3 normals_modelspace;

void main(void)
{
    //standard vertex positioning:
    gl_Position = mvpMatrix * vec4(vertex_modelspace,1);
    vertex_side = vertex_modelspace.x;

    //transformations for calculation of lighting:
    vec3 vertex_cameraspace=(vMatrix*mMatrix*vec4(vertex_modelspace,1)).xyz;

    //from vertex to camera (in cameraspace always at origin)
    EyeDirection_cameraspace = vec3(0,0,0) - vertex_cameraspace;

    //from vertex to light source
    LightDirection_cameraspace = LightPosition_cameraspace + EyeDirection_cameraspace;

    //Normals in camera space
    normals_cameraspace = (vMatrix * mMatrix * vec4(normals_modelspace,0)).xyz;
}
//! [0]
