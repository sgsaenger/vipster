#include "guidata.h"
#include "atom_model.h"
#include "bond_model.h"
#include <iostream>
#include <vector>
#ifndef __EMSCRIPTEN__
#include <QOpenGLContext>
#endif

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
{}

void GUI::GlobalData::initGL()
{
    if (initialized) return;
#ifndef __EMSCRIPTEN__
    initializeOpenGLFunctions();
    const auto fmt = QOpenGLContext::currentContext()->format();
    bool newEnough = false;
    if(QOpenGLContext::openGLModuleType() == QOpenGLContext::LibGL){
        if(fmt.version() >= qMakePair(3,3)){
            header = "# version 330\n";
            newEnough = true;
        }else{
            header = "# version 140\n";
        }
    }else{
        if(fmt.version() >= qMakePair(3,0)){
            header = "# version 300 es\nprecision highp float;\n";
            newEnough = true;
        }else{
            header = "# version 100 es\nprecision highp float;\n";
        }
    }
    folder = ":/shaders";
#else
    header = "# version 300 es\nprecision highp float;\n";
    folder = "";
#endif
    std::cout << "OpenGL Vendor: " << glGetString(GL_VENDOR) << std::endl;
    std::cout << "OpenGL Renderer: " << glGetString(GL_RENDERER) << std::endl;
    auto glVersionStr = reinterpret_cast<const char*>(glGetString(GL_VERSION));
    std::cout << "OpenGL Version: " << glVersionStr << std::endl;
    std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    if(!newEnough){
        std::cout << "Sorry, Vipster requires OpenGL 3.3 or OpenGL ES 3.0 or higher. Exiting." << std::endl;
        std::exit(1);
    }
    // skip ES information
    bool isES = strstr(glVersionStr, "OpenGL ES");
    if(isES){
        glVersionStr = glVersionStr+11;
    }
    char *glVersionEnd;
    int majorVersion = std::strtol(glVersionStr, &glVersionEnd, 10);
    int minorVersion = std::strtol(glVersionEnd+1, &glVersionEnd, 10);
    // generate buffers and upload data
    glGenBuffers(1, &sphere_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, sphere_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(atom_model),
                 static_cast<const void*>(&atom_model), GL_STATIC_DRAW);
    glGenBuffers(1, &cylinder_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, cylinder_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(bond_model),
                 static_cast<const void*>(&bond_model), GL_STATIC_DRAW);
    glGenBuffers(1, &cell_ibo);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cell_ibo);
    GLushort indices[24] = {0,1,0,2,0,3,1,4,1,5,2,4,2,6,3,5,3,6,4,7,5,7,6,7};
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), static_cast<void*>(indices), GL_STATIC_DRAW);
    initialized = true;
}

void GUI::Data::syncToGPU(void *context)
{
    #ifndef __EMSCRIPTEN__
    if(!wrap_initialized){
        initializeOpenGLFunctions();
        wrap_initialized = true;
    }
    #endif
    if(!initialized[context]){
        initGL(context);
        initialized[context] = true;
    }
    if(updated){
        updateGL();
        updated = false;
    }
}

GLuint GUI::Data::loadShader(const std::string &vert, const std::string &frag)
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
