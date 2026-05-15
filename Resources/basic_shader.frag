uniform sampler2D texture;
uniform float threshold;
uniform vec2 resolution;

void main() 
{
    vec2 uv = gl_TexCoord[0].xy;
    vec4 pixel = texture2D(texture, uv);
    
    //CALCULATE NORMALS for 3D effect
    //reminder our image colors are entirely from alpha so use alpha to represent height
    //check alpha values (height) of each pixel around current pixel

    vec2 offset = 1.0 / resolution;
    float left   = texture2D(texture, uv - vec2(offset.x, 0.0)).a;
    float right  = texture2D(texture, uv + vec2(offset.x, 0.0)).a;
    float top    = texture2D(texture, uv - vec2(0.0, offset.y)).a;
    float bottom = texture2D(texture, uv + vec2(0.0, offset.y)).a;

    //z param changes how sharp the lighting is. higher = flatter
    vec3 normal = normalize(vec3(right - left, bottom - top, 0.2)); 
        
    //LIGHTING
    vec3 lightDir = normalize(vec3(-1.0, -1.0, 1.5));
    float diffuse = max(dot(normal, lightDir), 0.0);
  
    vec3 viewDir = vec3(0.0, 0.0, 1.0);
    vec3 halfDir = normalize(lightDir + viewDir);
    float specular = pow(max(dot(normal, halfDir), 0.0), 16.0);

    //COLORING
    float depth = clamp((pixel.a - threshold) / (1.0 - threshold), 0.0, 1.0);
    vec3 deepColor = vec3(0.0, 0.25, 0.7);
    vec3 shallowColor = vec3(0.3, 0.8, 1.0);
    vec3 baseColor = mix(shallowColor, deepColor, depth);
 
    vec3 finalColor = baseColor * (0.5 + 0.5 * diffuse) + (specular * 0.4);

    //light blue halo around edge of particle
    float fresnel = pow(1.0 - max(dot(normal, viewDir), 0.0), 16.0);

    // Add a bright, glowing rim light to the edges
    vec3 rimColor = vec3(0.8, 0.95, 1.0); 
    finalColor += rimColor * fresnel * 0.6;
    //get rid of harsh edges
    float finalAlpha = smoothstep(threshold, threshold + 0.02, pixel.a);

    gl_FragColor = vec4(finalColor, finalAlpha);
}
