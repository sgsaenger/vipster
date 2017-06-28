#include "glwrapper.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <atom_model.h>
#include <bond_model.h>

std::string readShader(std::string filePath){
    std::string content;
    std::ifstream fileStream{filePath};

    if(!fileStream) throw std::invalid_argument{"Shader not found: "+filePath};

    content.assign(std::istreambuf_iterator<char>{fileStream},
                   std::istreambuf_iterator<char>{});
    return content;
}

GLuint loadShader(std::string header, std::string vertPath, std::string fragPath)
{
    GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    GLint gl_ok = GL_FALSE;

    std::string vertShaderStr = header+readShader(vertPath);
    const char *vertShaderSrc = vertShaderStr.c_str();
    glShaderSource(vertShader, 1, &vertShaderSrc, nullptr);
    glCompileShader(vertShader);
    glGetShaderiv(vertShader, GL_COMPILE_STATUS, &gl_ok);
    if(!gl_ok){
        std::cout << "Shader does not compile: " << vertPath << std::endl;
        GLint infoLen = 0;
        glGetShaderiv(vertShader, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog;
        infoLog.resize((infoLen > 1)?infoLen:1);
        glGetShaderInfoLog(vertShader, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Shader does not compile: "+vertPath};
    }

    std::string fragShaderStr = header+readShader(fragPath);
    const char *fragShaderSrc = fragShaderStr.c_str();
    glShaderSource(fragShader, 1, &fragShaderSrc, nullptr);
    glCompileShader(fragShader);
    glGetShaderiv(fragShader, GL_COMPILE_STATUS, &gl_ok);
    if(!gl_ok){
        std::cout << "Shader does not compile: " << fragPath << std::endl;
        GLint infoLen = 0;
        glGetShaderiv(fragShader, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog(infoLen > 1?infoLen:1);
        glGetShaderInfoLog(fragShader, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Shader does not compile: "+fragPath};
    }

    GLuint program = glCreateProgram();
    glAttachShader(program, vertShader);
    glAttachShader(program, fragShader);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &gl_ok);
    if(!gl_ok){
        std::cout << "Program does not link: " << vertPath << ", " << fragPath << std::endl;
        GLint infoLen = 0;
        std::vector<char> infoLog(infoLen > 1?infoLen:1);
        glGetProgramInfoLog(program, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Program does not link "+vertPath+" and "+fragPath};
    }

    glDetachShader(program, vertShader);
    glDeleteShader(vertShader);
    glDetachShader(program, fragShader);
    glDeleteShader(fragShader);
    return program;
}

void glMatScale(glMat& m, float f)
{
    for(int i=0;i<4;i++){
        m[i*4+0]*=f;
        m[i*4+1]*=f;
        m[i*4+2]*=f;
    }
}

GLWrapper::GLWrapper()
{
    glGenVertexArrays(1, &atom_vao);
    glBindVertexArray(atom_vao);
    glGenBuffers(1, &sphere_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    glBufferData(GL_ARRAY_BUFFER, atom_model_npoly*3*sizeof(float), (void*)&atom_model, GL_STATIC_DRAW);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,0);
    glEnableVertexAttribArray(0);
    glGenBuffers(1, &atom_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, atom_vbo);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,8*sizeof(float),0);
    glVertexAttribPointer(2,1,GL_FLOAT,GL_FALSE,8*sizeof(float),(const GLvoid*)(3*sizeof(float)));
    glVertexAttribPointer(3,4,GL_FLOAT,GL_FALSE,8*sizeof(float),(const GLvoid*)(4*sizeof(float)));
    glVertexAttribDivisor(1,1);
    glVertexAttribDivisor(2,1);
    glVertexAttribDivisor(3,1);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    glEnableVertexAttribArray(3);
}

GLWrapper::~GLWrapper()
{
    glDeleteBuffers(1, &sphere_vbo);
    glDeleteBuffers(1, &atom_vbo);
    glDeleteVertexArrays(1, &atom_vao);
}
