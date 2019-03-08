#include "guidata.h"
#include "atom_model.h"
#include "bond_model.h"
#include <iostream>
#include <vector>

using namespace Vipster;

#ifdef __EMSCRIPTEN__
#include <fstream>
std::string readShader(const std::string &filePath)
{
    std::string content;
    std::ifstream fileStream{filePath};

    if(!fileStream) throw std::invalid_argument{"Shader not found: "+filePath};

    content.assign(std::istreambuf_iterator<char>{fileStream},
                   std::istreambuf_iterator<char>{});
    return content;
}
#else
#include <QString>
#include <QFile>

std::string readShader(const std::string &filePath)
{
    QFile f(QString::fromStdString(filePath));
    f.open(QIODevice::ReadOnly);
    return f.readAll().toStdString();
}
#endif

GUI::GlobalData::GlobalData()
    : sphere_vbo{buffers[0]}, cylinder_vbo{buffers[1]},
      cell_ibo{buffers[2]}
{}

void GUI::GlobalData::initGL(const std::string& h, const std::string& f)
{
    header = h;
    folder = f;
#ifndef __EMSCRIPTEN__
    initializeOpenGLFunctions();
#endif
    // generate buffers and upload data
    glGenBuffers(3, buffers);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(atom_model),
                 static_cast<const void*>(&atom_model), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, cylinder_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(bond_model),
                 static_cast<const void*>(&bond_model), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_ibo);
    GLushort indices[24] = {0,1,0,2,0,3,1,4,1,5,2,4,2,6,3,5,3,6,4,7,5,7,6,7};
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), static_cast<void*>(indices), GL_STATIC_DRAW);
}

void GUI::Data::syncToGPU()
{
    #ifndef __EMSCRIPTEN__
    if(!wrap_initialized){
        initializeOpenGLFunctions();
        wrap_initialized = true;
    }
    #endif
    if(!initialized){
        initGL();
        initialized = true;
    }
    if(updated){
        updateGL();
        updated = false;
    }
}

GLuint GUI::Data::loadShader(std::string vert, std::string frag)
{
    GLuint vertShader = glCreateShader(GL_VERTEX_SHADER);
    GLuint fragShader = glCreateShader(GL_FRAGMENT_SHADER);
    GLint gl_ok = GL_FALSE;

    std::string vertShaderStr = global.header + readShader(global.folder + vert);
    const char *vertShaderSrc = vertShaderStr.c_str();
    glShaderSource(vertShader, 1, &vertShaderSrc, nullptr);
    glCompileShader(vertShader);
    glGetShaderiv(vertShader, GL_COMPILE_STATUS, &gl_ok);
    if(gl_ok == 0){
        std::cout << "Vertex-Shader does not compile" << std::endl;
        GLint infoLen = 0;
        glGetShaderiv(vertShader, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog;
        infoLog.resize((infoLen > 1)?static_cast<size_t>(infoLen):1);
        glGetShaderInfoLog(vertShader, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Vertex-Shader does not compile"};
    }

    std::string fragShaderStr = global.header + readShader(global.folder + frag);
    const char *fragShaderSrc = fragShaderStr.c_str();
    glShaderSource(fragShader, 1, &fragShaderSrc, nullptr);
    glCompileShader(fragShader);
    glGetShaderiv(fragShader, GL_COMPILE_STATUS, &gl_ok);
    if(gl_ok == 0){
        std::cout << "Fragment-Shader does not compile" << std::endl;
        GLint infoLen = 0;
        glGetShaderiv(fragShader, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog(infoLen > 1?static_cast<size_t>(infoLen):1);
        glGetShaderInfoLog(fragShader, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Fragment-Shader does not compile"};
    }

    GLuint program = glCreateProgram();
    glAttachShader(program, vertShader);
    glAttachShader(program, fragShader);
    glLinkProgram(program);
    glGetProgramiv(program, GL_LINK_STATUS, &gl_ok);
    if(gl_ok == 0){
        std::cout << "Program does not link" << std::endl;
        GLint infoLen = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &infoLen);
        std::vector<char> infoLog(infoLen > 1?static_cast<size_t>(infoLen):1);
        glGetProgramInfoLog(program, infoLen, nullptr, &infoLog[0]);
        std::cout << &infoLog[0] << std::endl;
        throw std::invalid_argument{"Program does not link"};
    }

    glDetachShader(program, vertShader);
    glDeleteShader(vertShader);
    glDetachShader(program, fragShader);
    glDeleteShader(fragShader);
    glUniformBlockBinding(program, 0, glGetUniformBlockIndex(program, "viewMat"));
    return program;
}

GUI::Data::Data(const GlobalData& glob)
    : global{glob}
{}

GUI::Data::Data(GUI::Data&& dat)
    : global{dat.global}
{
    std::swap(updated, dat.updated);
    std::swap(initialized, dat.initialized);
}
