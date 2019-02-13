in vec4 color_input;
out vec4 color_output;
flat in uint frag_discard;

void main(void)
{
    if(frag_discard!=uint(0)){
        discard;
    }else{
        color_output = color_input;
    }
}
