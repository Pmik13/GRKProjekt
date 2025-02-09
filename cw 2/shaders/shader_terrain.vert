#version 330 core

layout(location = 0) in vec3 inPosition;  // Vertex position (x, y, z)
layout(location = 1) in vec3 inNormal;    // Vertex normal
layout(location = 2) in vec3 inColor;     // Vertex color (RGB)

out vec3 FragPos;    // Position to pass to fragment shader
out vec3 FragNormal; // Normal to pass to fragment shader
out vec3 FragColor;  // Color to pass to fragment shader

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    // Calculate the position in world space
    FragPos = vec3(model * vec4(inPosition, 1.0));
    
    // Pass the normal to the fragment shader
    FragNormal = mat3(transpose(inverse(model))) * inNormal;
    
    // Pass the color to the fragment shader
    FragColor = inColor;

    // Apply transformations (model, view, projection matrices)
    gl_Position = projection * view * vec4(FragPos, 1.0);
}
