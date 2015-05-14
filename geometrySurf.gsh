#version 330

layout(points) in;
layout(triangle_strip,max_vertices=15) out;
//layout(points,max_vertices=15) out;

uniform vec3 shape;
uniform mat3 cellVec;
uniform mat4 vpMatrix;
uniform isamplerBuffer edgeOff;

in int vertices[1];
//in vec3 coord[1];
in int vID[1];

out vec4 pcol;

//Geometry:
// Bitmask for vertices:
//
//     3-------7
//    /|      /|
//   / |     / |
//  1--+----5  |
//  |  2----+--6
// z| /     | /y
//  |/      |/
//  0-------4
//      x
//
// Mask for edges:
//
//
//     +---9---+
//    /|      /|
//   4 8     5 10
//  +--+1---+  |
//  |  +--11+--+
//  0 /     2 /
//  |7      |6
//  +---3---+
//

int primPerCube[256] = int[256](
0, 1, 1, 2, 1, 2, 2, 3,
1, 2, 2, 3, 2, 3, 3, 2,
1, 2, 2, 3, 2, 3, 3, 4,
2, 3, 3, 4, 3, 4, 4, 3,
1, 2, 2, 3, 2, 3, 3, 4,
2, 3, 3, 4, 3, 4, 4, 3,
2, 3, 3, 2, 3, 4, 4, 3,
3, 4, 4, 3, 4, 3, 3, 2,
1, 2, 2, 3, 2, 3, 3, 4,
2, 3, 3, 4, 3, 4, 4, 3,
2, 3, 3, 4, 3, 2, 4, 3,
3, 4, 4, 3, 4, 3, 3, 2,
2, 3, 3, 4, 3, 4, 4, 3,
3, 4, 4, 3, 4, 3, 3, 2,
3, 4, 4, 3, 4, 3, 3, 2,
4, 3, 3, 2, 3, 2, 2, 1,
1, 2, 2, 3, 2, 3, 3, 4,
2, 3, 3, 4, 3, 4, 4, 3,
2, 3, 3, 4, 3, 4, 4, 3,
3, 4, 4, 3, 4, 3, 3, 2,
2, 3, 3, 4, 3, 4, 4, 3,
3, 4, 2, 3, 4, 3, 3, 2,
3, 4, 4, 3, 4, 3, 3, 2,
4, 3, 3, 2, 3, 2, 2, 1,
2, 3, 3, 4, 3, 4, 4, 3,
3, 4, 4, 3, 2, 3, 3, 2,
3, 4, 4, 3, 4, 3, 3, 2,
4, 3, 3, 2, 3, 2, 2, 1,
3, 4, 4, 3, 4, 3, 3, 2,
4, 3, 3, 2, 3, 2, 2, 1,
2, 3, 3, 2, 3, 2, 2, 1,
3, 2, 2, 1, 2, 1, 1, 0);

void main(void)
{
    pcol=vec4(0,0.5,0,0.5);

    vec3 coord;
    //coord.z = (0.5+mod(vID[0],shape.z))/shape.z;
    //coord.y = (0.5+mod(vID[0]/int(shape.z),shape.y))/shape.y;
    //coord.x = (0.5+mod(vID[0]/int(shape.z*shape.y),shape.x))/shape.x;
    coord.z = mod(vID[0],shape.z);
    coord.y = mod(vID[0]/int(shape.z),shape.y);
    coord.x = mod(vID[0]/int(shape.z*shape.y),shape.x);

    vec3 eoff[12];
    eoff[0]  = vec3( 0.0, 0.0, 0.5);
    eoff[1]  = vec3( 0.5, 0.0, 1.0);
    eoff[2]  = vec3( 1.0, 0.0, 0.5);
    eoff[3]  = vec3( 0.5, 0.0, 0.0);
    eoff[4]  = vec3( 0.0, 0.5, 1.0);
    eoff[5]  = vec3( 1.0, 0.5, 1.0);
    eoff[6]  = vec3( 1.0, 0.5, 0.0);
    eoff[7]  = vec3( 0.0, 0.5, 0.0);
    eoff[8]  = vec3( 0.0, 1.0, 0.5);
    eoff[9]  = vec3( 0.5, 1.0, 1.0);
    eoff[10] = vec3( 1.0, 1.0, 0.5);
    eoff[11] = vec3( 0.5, 1.0, 0.0);
    //eoff[0]  = vec3(-0.5,-0.5, 0.0);
    //eoff[1]  = vec3( 0.0,-0.5, 0.5);
    //eoff[2]  = vec3( 0.5,-0.5, 0.0);
    //eoff[3]  = vec3( 0.0,-0.5,-0.5);
    //eoff[4]  = vec3(-0.5, 0.0, 0.5);
    //eoff[5]  = vec3( 0.5, 0.0, 0.5);
    //eoff[6]  = vec3( 0.5, 0.0,-0.5);
    //eoff[7]  = vec3(-0.5, 0.0,-0.5);
    //eoff[8]  = vec3(-0.5, 0.5, 0.0);
    //eoff[9]  = vec3( 0.0, 0.5, 0.5);
    //eoff[10] = vec3( 0.5, 0.5, 0.0);
    //eoff[11] = vec3( 0.0, 0.5,-0.5);

    for (int i=0;i<primPerCube[vertices[0]];i++)
    {
        int j;
        j = (12*vertices[0])+(i*3);
        int e1,e2,e3;
        e1 = texelFetch(edgeOff,j).r;
        e2 = texelFetch(edgeOff,j+1).r;
        e3 = texelFetch(edgeOff,j+2).r;
        vec3 c1,c2,c3;
        c1 = (coord+eoff[e1])/shape;
        c2 = (coord+eoff[e2])/shape;
        c3 = (coord+eoff[e3])/shape;
        //c1 = coord+edgeOff[j];
        //c2 = coord+edgeOff[j+1];
        //c3 = coord+edgeOff[j+2];

        gl_Position = vpMatrix * (vec4(c1*cellVec,1)+gl_in[0].gl_Position);
        EmitVertex();
        gl_Position = vpMatrix * (vec4(c2*cellVec,1)+gl_in[0].gl_Position);
        EmitVertex();
        gl_Position = vpMatrix * (vec4(c3*cellVec,1)+gl_in[0].gl_Position);
        EmitVertex();
        EndPrimitive();
    }
}
