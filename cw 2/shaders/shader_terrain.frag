#version 330 core
out vec4 FragColor;

in vec3 FragPos;  // Receive the position from the vertex shader

void main() {
    // Set the color to blue
    FragColor = vec4(0.0, 0.0, 1.0, 1.0); // Blue color
}
