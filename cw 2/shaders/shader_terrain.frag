#version 330 core

in vec3 vertexColor;  // Input from vertex shader
out vec4 FragColor;   // Final color output to the screen

void main() {
    FragColor = vec4(vertexColor, 1.0);  // Set the color of the fragment (pixel)
}
