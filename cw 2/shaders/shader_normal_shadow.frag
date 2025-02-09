#version 430 core

float AMBIENT = 0.1;

uniform vec3 lightPos;
uniform sampler2D colorTexture;
uniform sampler2D normalSampler;
uniform sampler2D depthMap;
uniform vec3 lightDir;

in vec2 vecTex;
in vec3 viewDirTS;
in vec3 lightDirTS;
in vec4 shadowCoord; // Pozycja w przestrzeni œwiat³a
in vec3 normal;

out vec4 outColor;

// Funkcja sprawdzaj¹ca cieñ
float calculateShadow(vec4 shadowCoord)
{
    vec3 projCoords = shadowCoord.xyz / shadowCoord.w;
    projCoords = projCoords * 0.5 + 0.5; // Przekszta³cenie do zakresu [0,1]

    float closestDepth = texture(depthMap, projCoords.xy).r; // G³êbokoœæ z mapy cieni
    float currentDepth = projCoords.z;

    float bias = max(0.0025 * (1.0 - dot(normal, lightDir)), 0.0005);
    float shadow = currentDepth - bias > closestDepth ? 0.3 : 1.0; // 0.3 oznacza cieñ

    return shadow;
}

void main()
{
    vec3 lightDir = normalize(lightDirTS);
    vec3 viewDir = normalize(viewDirTS);
    
    vec3 textureColor = texture(colorTexture, vecTex).rgb;
    vec3 normal = normalize(texture(normalSampler, vecTex).xyz * 2.0 - 1.0);

    float diffuse = max(dot(normal, lightDir), 0.0);
    float shadow = calculateShadow(shadowCoord);

    vec3 finalColor = textureColor * (AMBIENT + diffuse) * shadow;
    
    outColor = vec4(finalColor, 1.0);
}
