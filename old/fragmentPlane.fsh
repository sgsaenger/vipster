#version 330

out vec4 planeColor;
in vec2 UV;
uniform sampler2D texSampler;

void main(void)
{
    float v = texture(texSampler,UV).r;
    float red = exp(-6.25*pow(v-0.9,2));
    float green = exp(-6.25*pow(v-0.5,2));
    float blue = exp(-6.25*pow(v-0.1,2));
    planeColor = vec4(red,green,blue,0.5);
}
