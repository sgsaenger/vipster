#include "glwrapper.h"
#include <fstream>
#include <iostream>
#include <vector>

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

    glDeleteShader(vertShader);
    glDeleteShader(fragShader);
    return program;
}
