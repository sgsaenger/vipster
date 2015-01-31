#version 330

//! [0]
layout(location = 0) in vec3 position_modelspace;
layout(location=1) in float scale_modelspace;
in vec3 vertex_modelspace;
in vec3 normals_modelspace;

uniform mat4 mvpMatrix;
uniform mat4 mvMatrix;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;


void main(void)
{
    //standard vertex positioning:
    // CHANGE TO: mvpMatrix * vec4(vertex_modelspace*scale + pos_shift,1);
    gl_Position = mvpMatrix * vec4(vertex_modelspace*scale_modelspace+position_modelspace,1);
    //gl_Position=mvpMatrix*vec4(vertex_modelspace+position_modelspace,1);

    //transformations for calculation of lighting:
    vec3 vertex_cameraspace=(mvMatrix*vec4(vertex_modelspace,1)).xyz;

    //from vertex to camera (in cameraspace always at origin)
    EyeDirection_cameraspace = vec3(0,0,0) - vertex_cameraspace;

    //from vertex to light source
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;

    //Normals in camera space
    normals_cameraspace = (mvMatrix * vec4(normals_modelspace,0)).xyz;
}
//! [0]
