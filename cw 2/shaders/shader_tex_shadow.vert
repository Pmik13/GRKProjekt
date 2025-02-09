#version 430 core

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexNormal;
layout(location = 2) in vec2 vertexTexCoord;

uniform mat4 transformation;
uniform mat4 modelMatrix;
uniform mat4 LightVP; // Macierz œwiat³a

out vec3 vecNormal;
out vec3 worldPos;
out vec2 vecTex;
out vec4 shadowCoord; // Pozycja w przestrzeni œwiat³a

void main()
{
    worldPos = (modelMatrix * vec4(vertexPosition, 1)).xyz;
    vecNormal = normalize(mat3(modelMatrix) * vertexNormal);
    vecTex = vec2(vertexTexCoord.x, 1.0 - vertexTexCoord.y);
    
    shadowCoord = LightVP * vec4(worldPos, 1.0); // Transformacja do przestrzeni œwiat³a
    gl_Position = transformation * vec4(vertexPosition, 1.0);
}