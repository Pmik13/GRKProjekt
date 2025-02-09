#version 430 core

float AMBIENT = 0.4;

uniform vec3 color;
uniform vec3 lightPos;
uniform sampler2D colorTexture;
uniform sampler2D depthMap; // Mapa cieni
uniform vec3 lightDir;

in vec3 vecNormal;
in vec3 worldPos;
in vec2 vecTex;
in vec4 shadowCoord; // Pozycja w przestrzeni œwiat³a

out vec4 outColor;

// Funkcja sprawdzaj¹ca, czy piksel jest w cieniu
float calculateShadow(vec4 shadowCoord)
{
    vec3 projCoords = shadowCoord.xyz / shadowCoord.w;
    projCoords = projCoords * 0.5 + 0.5; // Przekszta³cenie do zakresu [0,1]

    float closestDepth = texture(depthMap, projCoords.xy).r; // G³êbokoœæ z mapy cieni
    float currentDepth = projCoords.z;

    float bias = max(0.0025 * (1.0 - dot(vecNormal, lightDir)), 0.0005);
    float shadow = currentDepth - bias > closestDepth ? 0.4 : 1.0; // 0.4 to czêœciowy cieñ

    return shadow;
}

void main()
{
    vec3 lightDir = normalize(lightPos - worldPos);
    vec3 normal = normalize(vecNormal);
    vec3 textureColor = texture(colorTexture, vecTex).rgb;

    float diffuse = max(dot(normal, lightDir), 0.0);
    
    float shadow = calculateShadow(shadowCoord);
    
    vec3 finalColor = textureColor * min(1.0, AMBIENT + diffuse) * shadow;

    outColor = vec4(finalColor, 1.0);
}
