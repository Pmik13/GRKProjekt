#version 430 core

float AMBIENT = 0.1;

uniform vec3 color;
uniform vec3 lightPos;
uniform sampler2D colorTexture;
uniform sampler2D normalSampler;

in vec2 vecTex;
in vec3 viewDirTS;
in vec3 lightDirTS;

out vec4 outColor;
void main()
{
	vec3 lightDir = normalize(lightDirTS);
	vec3 viewDir = normalize(viewDirTS);

	vec3 textureColor = texture2D(colorTexture, vecTex).xyz;
	vec3 normal = normalize((texture2D(normalSampler, vecTex).xyz)*2.0 - 1.0);

	float diffuse=max(0,dot(normal,lightDir));
	outColor = vec4(textureColor*min(1,AMBIENT+diffuse), 1.0);
}