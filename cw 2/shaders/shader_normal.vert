#version 430 core

layout(location = 0) in vec3 vertexPosition;
layout(location = 1) in vec3 vertexNormal;
layout(location = 2) in vec2 vertexTexCoord;
layout(location = 3) in vec3 vertexTangent;
layout(location = 4) in vec3 vertexBitangent;

uniform mat4 transformation;
uniform mat4 modelMatrix;
uniform mat4 LightVP;
uniform vec3 lightPos;
uniform vec3 cameraPos;

out vec2 vecTex;
out vec3 viewDirTS;
out vec3 lightDirTS;
out vec4 shadowCoord; // Pozycja w przestrzeni œwiat³a
out vec3 normal;

void main()
{
    vec3 worldPos = (modelMatrix * vec4(vertexPosition, 1)).xyz;
    vec3 normal = normalize((modelMatrix * vec4(vertexNormal, 0.0)).xyz);
    vec3 tangent = normalize((modelMatrix * vec4(vertexTangent, 0.0)).xyz);
    vec3 bitangent = normalize((modelMatrix * vec4(vertexBitangent, 0.0)).xyz);

    mat3 TBN = transpose(mat3(tangent, bitangent, normal));

    vec3 viewDir = normalize(cameraPos - worldPos);
    vec3 lightDir = normalize(lightPos - worldPos);

    viewDirTS = normalize(TBN * viewDir);
    lightDirTS = normalize(TBN * lightDir);

    vecTex = vec2(vertexTexCoord.x, 1.0 - vertexTexCoord.y);
    shadowCoord = LightVP * vec4(worldPos, 1.0); // Przekszta³cenie do przestrzeni œwiat³a

    gl_Position = transformation * vec4(vertexPosition, 1.0);
}