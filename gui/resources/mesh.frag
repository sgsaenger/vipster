in vec3 normals_cameraspace;
in vec3 EyeDirection_cameraspace;
in vec3 LightDirection_cameraspace;
in vec2 UV;

uniform sampler2D tex;

out vec4 color_output;

void main(void)
{
    vec4 color_input = texture(tex, UV);
    //hardcoded lighting parameters.
    //think about putting them as uniforms.
    vec3 LightColor = vec3(1,1,1);
    float LightPower = 50.0f;
    float distance = 10.;

    //Material properties:
    vec3 MaterialAmbientColor = vec3(0.5,0.5,0.5)*color_input.xyz;
    vec3 MaterialSpecularColor = vec3(0.3,0.3,0.3);

    //diffuse lighting:
    vec3 n = normalize(normals_cameraspace);
    vec3 l = normalize(LightDirection_cameraspace);
    //angle of incidence
    float cosTheta = clamp(dot(n,l),0.,1.);

    //specular lighting:
    vec3 E = normalize(EyeDirection_cameraspace);
    //direct reflection vector:
    vec3 R = reflect(-l,n);
    //angular difference from R:
    float cosAlpha = clamp(dot(E,R),0.,1.);

    //final color
    vec3 fragTemp = MaterialAmbientColor +
                MaterialSpecularColor*LightColor*LightPower * pow(cosAlpha,10.)/(distance*distance)+
                color_input.xyz*LightColor*LightPower*cosTheta/(distance*distance);
    color_output = vec4(fragTemp,color_input.a);
}
