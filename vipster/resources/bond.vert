layout(location = 0) in vec3 vertex;
layout(location = 1) in mat3 mMatrix;
layout(location = 4) in vec3 position;
layout(location = 5) in ivec4 pbc_crit;
layout(location = 6) in vec4 color1;
layout(location = 7) in vec4 color2;

layout(std140, row_major) uniform viewMat{
    mat4 vpMatrix;
    mat4 rMatrix;
};
uniform mat3 pos_scale;
uniform vec3 offset;
uniform ivec3 pbc_cell;
uniform ivec3 mult;

out vec3 normals_cameraspace;
out vec3 EyeDirection_cameraspace;
out vec3 LightDirection_cameraspace;
out float vertex_side;
out vec4 s1Cpass;
out vec4 s2Cpass;
flat out uint render;

void main(void)
{
    if(pbc_crit.w != 0){
        render = uint(0);
    }else{
        ivec3 test = ivec3(mult) - ivec3(pbc_crit.xyz);
        render = uint((test.x > pbc_cell.x) && (test.y > pbc_cell.y) && (test.z > pbc_cell.z));
    }
    if(render!=uint(0)){
        //standard vertex positioning:
        gl_Position = vpMatrix * vec4(mMatrix * vertex + position * pos_scale + offset, 1);
        //pass coordinate and colors to fragment shader
        vertex_side = vertex.x;
        s1Cpass = color1;
        s2Cpass = color2;

        //transformations for calculation of lighting:
        vec3 vertex_cameraspace=(rMatrix*vec4(mMatrix*vertex,0)).xyz;

        //from vertex to camera (in cameraspace always at origin)
        EyeDirection_cameraspace = vec3(0,0,25) - vertex_cameraspace;

        //from vertex to light source
        LightDirection_cameraspace = vec3(10,10,10) + EyeDirection_cameraspace;

        //Normals in camera space
        normals_cameraspace = (rMatrix * vec4(mMatrix * vec3(0,vertex.yz),0)).xyz;
    }
}
