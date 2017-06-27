#ifndef GLWRAPPER_H
#define GLWRAPPER_H

#include <string>
#include <GLES3/gl3.h>
#include <vector>
#include <array>

GLuint loadShader(std::string header, std::string vertPath, std::string fragPath);

class GLWrapper{
    public:
        GLWrapper();
        ~GLWrapper();
        GLuint atom_program, bond_program, cell_program;
        GLuint atom_vao, bond_vao, cell_vao;
        GLuint atom_vbo, bond_vbo, cell_vbo;
        GLuint sphere_vbo, torus_vbo;
        std::vector<std::array<float,8>> atom_buffer;
};

#endif // GLWRAPPER_H
