out vec4 fragColor;

in vec3 normals_cameraspace;
in vec3 EyeDirection_cameraspace;
in vec3 LightDirection_cameraspace;
in float vertex_side;
in vec4 s1Cpass;
in vec4 s2Cpass;
flat in uint render;

void main(void)
{
    if(render==uint(0)){
        discard;
    }else{
        //determine color:
        vec4 MaterialDiffuseColor = vec4(0,0,0,1);
        if (vertex_side > 0.0)
        {
            MaterialDiffuseColor = s1Cpass;
        }
        else
        {
            MaterialDiffuseColor = s2Cpass;
        }
        //hardcoded lighting parameters.
        //think about putting them as uniforms.
        vec3 LightColor = vec3(1,1,1);
        float LightPower = 50.0f;
        float distance = 10.;

        //Material properties:
        vec3 MaterialAmbientColor = vec3(0.5,0.5,0.5)*MaterialDiffuseColor.rgb;
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
                    MaterialDiffuseColor.rgb*LightColor*LightPower*cosTheta/(distance*distance);
        fragColor = vec4(fragTemp,MaterialDiffuseColor.a);
    }
}
