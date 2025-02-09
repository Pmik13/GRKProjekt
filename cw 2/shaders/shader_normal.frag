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

// Function to calculate shadow
float calculateShadow(vec4 shadowCoord)
{
    // Normalize shadow coordinates and transform to texture space (range [0, 1])
    vec3 projCoords = shadowCoord.xyz / shadowCoord.w;
    projCoords = projCoords * 0.5 + 0.5; // Transform to [0, 1]

    // If the projected Z value is outside [0, 1], there’s no shadow
    if (projCoords.z > 1.0)
        return 1.0;

    // Sample the depth map at the projected coordinates
    float closestDepth = texture(depthMap, projCoords.xy).r;
    float currentDepth = projCoords.z;

    // Compute a bias to avoid shadow acne
    float bias = max(0.0025 * (1.0 - dot(normal, lightDir)), 0.0005);

    // Apply Percentage Closer Filtering (PCF) for smoother shadow edges
    float shadow = 0.0;
    vec2 texelSize = 1.0 / textureSize(depthMap, 0);  // Size of a single texel in the shadow map

    // Loop to sample 9 neighboring texels around the current texel
    for (int x = -1; x <= 1; x++) {
        for (int y = -1; y <= 1; y++) {
            // Get the depth of the neighboring texels
            float pcfDepth = texture(depthMap, projCoords.xy + vec2(x, y) * texelSize).r;
            
            // Compare the current depth with the depth from the shadow map
            shadow += currentDepth - bias > pcfDepth ? 0.3 : 1.0;
        }
    }

    shadow /= 9.0f; // Average of the 9 samples to get the final shadow value

    return shadow;
}

void main()
{
    vec3 lightDirNorm = normalize(lightDirTS);  // Correct the naming to avoid confusion
    vec3 viewDir = normalize(viewDirTS);
    
    vec3 textureColor = texture(colorTexture, vecTex).rgb;
    vec3 normal = normalize(texture(normalSampler, vecTex).xyz * 2.0 - 1.0);

    float diffuse = max(dot(normal, lightDirNorm), 0.0);  // Use lightDirNorm
    float shadow = calculateShadow(shadowCoord);

    vec3 finalColor = textureColor * (AMBIENT + diffuse) * shadow;
    
    outColor = vec4(finalColor, 1.0);
}
