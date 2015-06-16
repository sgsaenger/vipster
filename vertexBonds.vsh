in vec3 vertex_modelspace;
#if __VERSION__ < 330
uniform mat4 mMatrix;
uniform vec4 s1Color;
uniform vec4 s2Color;
#else
in mat4 mMatrix;
in vec4 s1Color;
in vec4 s2Color;
#endif

uniform mat4 vpMatrix;
uniform mat4 rMatrix;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out float vertex_side;
out vec4 s1Cpass;
out vec4 s2Cpass;

void main(void)
{
    //standard vertex positioning:
    gl_Position = vpMatrix * mMatrix * vec4(vertex_modelspace,1);
    //pass coordinate and colors to fragment shader
    vertex_side = vertex_modelspace.x;
    s1Cpass = s1Color;
    s2Cpass = s2Color;

    //transformations for calculation of lighting:
    vec3 vertex_cameraspace=(rMatrix*mMatrix*vec4(vertex_modelspace,0)).xyz;

    //from vertex to camera (in cameraspace always at origin)
    EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;

    //from vertex to light source
    LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;

    //Normals in camera space
    normals_cameraspace = (rMatrix*mMatrix * vec4(0,vertex_modelspace.yz,0)).xyz;
}
