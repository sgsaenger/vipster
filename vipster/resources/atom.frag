out vec4 fragColor;

in float radius;
in vec2 center;
in vec4 MaterialDiffuseColor;
flat in uint frag_discard;

layout(std140, row_major) uniform viewMat{
    mat4 rvpMatrix;
    mat4 pMatrix;
    mat4 vMatrix;
    mat4 rMatrix;
};
uniform vec4 viewport;

void main(void)
{
    if(frag_discard!=uint(0)){
        discard;
    }else{
        // map back from screen space to clip space
        vec2 clipPos = (((2*gl_FragCoord.xy) - (2*viewport.xy)) / (viewport.zw)) - 1;
        // map back from clip space to view space
        vec2 viewPos = vec2((clipPos.x - pMatrix[0][3])/pMatrix[0][0],
                            (clipPos.y - pMatrix[1][3])/pMatrix[1][1]);
        // TODO: precalculate in vertex-shader?
        vec2 viewCenter = vec2((center.x - pMatrix[0][3])/pMatrix[0][0],
                               (center.y - pMatrix[1][3])/pMatrix[1][1]);

        vec2 diff = viewPos - viewCenter;
        float diff2 = dot(diff, diff);
        float radius2 = radius * radius;
        if (diff2 > radius2){
            // TODO: manual antialiasing
            discard;
        }else{
            float dr = sqrt(radius2 - diff2);
            // calculate lighting
//            vec3 MaterialAmbientColor = 0.5 * MaterialDiffuseColor.xyz;
//            vec3 MaterialSpecularColor = vec3(0.3, 0.3, 0.3);
//            vec3 fragTemp = MaterialAmbientColor;
            vec3 l = normalize(vec3(10, 10, 10));
            vec3 n = vec3(diff, dr);
            float intensity = .2 + max(dot(l, normalize(n)), 0.0);
            vec3 fragTemp = intensity * MaterialDiffuseColor.xyz;
            fragColor = vec4(fragTemp, MaterialDiffuseColor.a);
//            fragColor = MaterialDiffuseColor;
            // calculate z-position
            gl_FragDepth = gl_FragCoord.z + dr*pMatrix[2][2] + pMatrix[2][3];
        }
    }
}
