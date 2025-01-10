#include "glew.h"
#include <GLFW/glfw3.h>
#include "glm.hpp"
#include "ext.hpp"
#include <iostream>
#include <cmath>

#include "Shader_Loader.h"
#include "Render_Utils.h"
#include "Camera.h"

#include "Box.cpp"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <string>



GLuint program;
Core::Shader_Loader shaderLoader;
GLuint VAO, VBO;
Core::RenderContext shipContext;
Core::RenderContext sphereContext;
float Rotation1 = 0.0f;
float Rotation2 = 0.0f;
float Rotation3 = 0.0f;
float Rotation4 = 0.0f;
float Rotation5 = 0.0f;
float cameraAngle = 0;
glm::vec3 cameraPos = glm::vec3(-5, 0, 0);
glm::vec3 cameraDir;


glm::mat4 createCameraMatrix()
{
	cameraDir = glm::vec3(cosf(cameraAngle), 0.0f, sinf(cameraAngle));
	glm::vec3 up = glm::vec3(0, 1, 0);

	return Core::createViewMatrix(cameraPos, cameraDir, up);
}
void renderScene(GLFWwindow* window)
{
	glClearColor(0.0f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(program);

	// Kamera i perspektywa
	glm::mat4 camera = createCameraMatrix();
	glm::mat4 perspective = Core::createPerspectiveMatrix();

	glm::mat4 objectModel = glm::mat4(1.0f);

	// Przesunięcie obiektu do początku układu współrzędnych
	objectModel = glm::translate(objectModel, glm::vec3(0.0f, 0.0f, Rotation2 / 20));

	objectModel = glm::translate(objectModel, glm::vec3(0.0f, -5.4f, 2.5f));

	// Rotacja obiektu wokół własnej osi
	objectModel = glm::rotate(objectModel, glm::radians(Rotation1), glm::vec3(0.0f, 1.0f, 0.0f));

	// Powrót obiektu do pierwotnej pozycji
	objectModel = glm::translate(objectModel, glm::vec3(0.0f, 5.4f, -2.5f));

	// Translacja w kierunku globalnym



	glm::mat4 objectModelY = objectModel;
	objectModelY = glm::translate(objectModelY, glm::vec3(0.0f, Rotation3/50, 0.0f));

	glm::mat4 objectModelX = objectModelY;
	objectModelX = glm::translate(objectModelX, glm::vec3(0.0f, 1.1f, 1.5f));
	objectModelX = glm::rotate(objectModelX, glm::radians(Rotation5), glm::vec3(1.0f, 0.0f, 0.0f));
	objectModelX = glm::translate(objectModelX, glm::vec3(0.0f, -1.1f, -1.5f));


	glm::mat4 objectModelC = objectModelX;
	objectModelC  = glm::translate(objectModelC, glm::vec3(0.0f, +3.0f, 0.6f));
	objectModelC = glm::rotate(objectModelC, glm::radians(-Rotation4), glm::vec3(1.0f, 0.0f, 0.0f));
	objectModelC = glm::translate(objectModelC, glm::vec3(0.0f, -3.0f, -0.6f));




	// Renderowanie pierwszego prostokąta
	glm::mat4 model1 = glm::translate(objectModelC, glm::vec3(0.0f, 3.0f, 0.0f)); // Prostokąt w centrum
	model1 = glm::rotate(model1, glm::radians(-0.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	model1 = glm::scale(model1, glm::vec3(1.0f, 0.5f, 1.5f));
	model1 = glm::translate(model1, glm::vec3(0.0f, 0.0f, -0.3f));

	glm::mat4 transformation1 = perspective * camera * model1;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, glm::value_ptr(transformation1));
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	glm::mat4 model2 = glm::translate(objectModelX, glm::vec3(0.0f, 2.2f, 0.8f)); // Prostokąt w centrum
	model2 = glm::rotate(model2, glm::radians(60.0f), glm::vec3(1.0f, 0.0f, 0.0f));

	model2 = glm::scale(model2, glm::vec3(1.0f, 0.5f, 2.0f));

	glm::mat4 transformation2 = perspective * camera * model2;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, glm::value_ptr(transformation2));
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);

	glm::mat4 model3 = glm::translate(objectModelY, glm::vec3(0.0f, 0.2f, 1.5f)); // Prostokąt w centrum
	model3 = glm::rotate(model3, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));

	model3 = glm::scale(model3, glm::vec3(1.0f, 0.5f, 2.5f));

	glm::mat4 transformation3 = perspective * camera * model3;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, glm::value_ptr(transformation3));
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);


	glm::mat4 model4 = glm::translate(objectModelY, glm::vec3(0.0f, -1.0f, 2.0f)); // Prostokąt w centrum
	model4 = glm::rotate(model4, glm::radians(0.0f), glm::vec3(1.0f, 0.0f, 0.0f));

	model4 = glm::scale(model4, glm::vec3(1.0f, 0.5f, 1.0f));

	glm::mat4 transformation4 = perspective * camera * model4;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, glm::value_ptr(transformation4));
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);


	glm::mat4 model5 = glm::translate(objectModel, glm::vec3(0.0f, -1.0f, 2.5f)); // Prostokąt w centrum
	model5 = glm::rotate(model5, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	model5 = glm::scale(model5, glm::vec3(1.0f, 0.5f, 2.5f));


	glm::mat4 transformation5 = perspective * camera * model5;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, glm::value_ptr(transformation5));
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);


	glm::mat4 model6 = glm::translate(objectModel, glm::vec3(0.0f, -2.7f, 2.5f)); // Prostokąt w centrum
	model6 = glm::rotate(model6, glm::radians(90.0f), glm::vec3(1.0f, 0.0f, 0.0f));
	model6 = glm::scale(model6, glm::vec3(1.0f, 2.0f, 1.0f));


	glm::mat4 transformation6 = perspective * camera * model6;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, glm::value_ptr(transformation6));
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);


	glUseProgram(0);

	glfwSwapBuffers(window);
}
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

void init(GLFWwindow* window)
{
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	glEnable(GL_DEPTH_TEST);
	program = shaderLoader.CreateProgram("shaders/shader_2_2.vert", "shaders/shader_2_2.frag");

	float cubeVertices[] = {
		// Pozycja           // Kolor
		// Przednia ściana
		-0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,  // Czerwony
		 0.5f, -0.5f,  0.5f,  0.0f, 1.0f, 0.0f,  // Zielony
		-0.5f,  0.5f,  0.5f,  0.0f, 0.0f, 1.0f,  // Niebieski
		 0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,  // Żółty

		 // Tylna ściana
		 -0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 1.0f,  // Magenta
		  0.5f, -0.5f, -0.5f,  0.0f, 1.0f, 1.0f,  // Cyan
		 -0.5f,  0.5f, -0.5f,  0.5f, 0.5f, 0.5f,  // Szary
		  0.5f,  0.5f, -0.5f,  0.0f, 0.0f, 0.0f,  // Czarny
	};

	unsigned int cubeIndices[] = {
		// Przednia ściana
		0, 1, 2,
		1, 3, 2,

		// Tylna ściana
		4, 5, 6,
		5, 7, 6,

		// Lewa ściana
		0, 2, 4,
		2, 6, 4,

		// Prawa ściana
		1, 5, 3,
		5, 7, 3,

		// Górna ściana
		2, 3, 6,
		3, 7, 6,

		// Dolna ściana
		0, 1, 4,
		1, 5, 4,
	};

	// Generowanie i konfiguracja VAO/VBO
	unsigned int EBO;
	glGenVertexArrays(1, &VAO);
	glBindVertexArray(VAO);

	// Wczytanie danych wierzchołków
	glGenBuffers(1, &VBO);
	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, sizeof(cubeVertices), cubeVertices, GL_STATIC_DRAW);

	// Wczytanie danych indeksów
	glGenBuffers(1, &EBO);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(cubeIndices), cubeIndices, GL_STATIC_DRAW);

	// Konfiguracja atrybutów wierzchołków
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0); // Pozycje
	glEnableVertexAttribArray(0);

	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float))); // Kolory
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);

}

void shutdown(GLFWwindow* window)
{
	shaderLoader.DeleteProgram(program);
}

//obsluga wejscia
void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}


	float angleSpeed = 0.001f;
	float angleSpeed3 = 0.1f;
	float angleSpeed2 = 0.1f;
	float moveSpeed = 0.001f;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
		cameraAngle -= angleSpeed; 
	}
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
		cameraAngle += angleSpeed; 
	}
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		cameraPos += cameraDir * moveSpeed;
	}
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
		cameraPos -= cameraDir * moveSpeed;
	}
	if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
		Rotation1 -= angleSpeed2;
	}
	if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
		Rotation2 -= angleSpeed2;
	}
	if (glfwGetKey(window, GLFW_KEY_U) == GLFW_PRESS) {
		Rotation2 += angleSpeed2;
	}
	if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
		Rotation3 -= angleSpeed3;
	}
	if (glfwGetKey(window, GLFW_KEY_T) == GLFW_PRESS) {
		Rotation3 += angleSpeed3;
	}
	if (glfwGetKey(window, GLFW_KEY_C) == GLFW_PRESS) {
		Rotation4 += angleSpeed3;
	}
	if (glfwGetKey(window, GLFW_KEY_X) == GLFW_PRESS) {
		Rotation5 += angleSpeed3;
	}
}

// funkcja jest glowna petla
void renderLoop(GLFWwindow* window) {
	while (!glfwWindowShouldClose(window))
	{
		processInput(window);

		renderScene(window);
		glfwPollEvents();
	}
}
//}