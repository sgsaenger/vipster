#ifndef GLWRAPPER_H
#define GLWRAPPER_H

#include <string>
#include <GLES3/gl3.h>

GLuint loadShader(std::string header, std::string vertPath, std::string fragPath);

#endif // GLWRAPPER_H
