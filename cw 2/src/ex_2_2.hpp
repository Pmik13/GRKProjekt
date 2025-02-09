#include "glew.h"
#include <GLFW/glfw3.h>
#include "glm.hpp"
#include "ext.hpp"
#include <iostream>
#include <cmath>
#include <chrono>

#include "Shader_Loader.h"
#include "Render_Utils.h"
#include "stb_image.h"
#include "Camera.h"
#include "Texture.h"
#include "FastNoiseLite.h"

#include "Box.cpp"
#include <assimp/Importer.hpp>
#include <assimp/scene.h>
#include <assimp/postprocess.h>
#include <string>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"


GLuint program;
GLuint programTexShadow;
GLuint programNormalShadow;
GLuint programTex;
GLuint programNormal;
GLuint programTerrain;
Core::Shader_Loader shaderLoader;
GLuint shadowShaderProgram;
GLuint VAO, VBO;
GLuint VAOTer, VBOTer, EBOTer;
GLuint depthMapFBO;
GLuint depthMap;
GLuint skyboxTexture;
GLuint programSkybox;
Core::RenderContext skyboxContext;
Core::RenderContext shipContext;
Core::RenderContext buildingContext;
Core::RenderContext coneContext;
Core::RenderContext sphereContext;
float Rotation1 = 0.0f;
float Rotation2 = 0.0f;
float Rotation3 = 0.0f;
float cameraAngle = 0;
float maxSpeed = 1.0f;
bool isPopupOpen = false;
bool isButtonClicked = false;
std::chrono::steady_clock::time_point popupOpenTime;

glm::vec3 obstacleboxsize = glm::vec3(0.05f, 0.05f, 0.05f);
glm::vec3 buildingboxsize = glm::vec3(0.5f, 0.5f, 0.5f);
glm::vec3 boxsize = glm::vec3(0.2f, 0.2f, 0.2f);
glm::vec3 cameraPos = glm::vec3(-4.f, 0.0f, 0.0f);
glm::vec3 cameraDir = glm::vec3(1.f, 0.f, 0.f);;
glm::vec3 spaceshipPos = glm::vec3(-4.f, 0, 0);
glm::vec3 spaceshipDir = glm::vec3(1.f, 0.f, 0.f);

float near_plane = 0.1;
float far_plane = 40.;
glm::vec3 lightPos = glm::vec3(-5.0f, 30.0f, 5.0f);  // Światło nad sceną
glm::vec3 lightDir = glm::vec3(-1.0f, -1.5f, -1.0f);  // Kierunek w dół pod kątem

glm::mat4 lightVP = glm::ortho(-100.f, 100.f, -30.f, 30.f, near_plane, far_plane) *
glm::lookAt(lightPos, lightPos + lightDir, glm::vec3(0, 1, 0));

bool shadowMappingEnabled = true;
bool normalMappingEnabled = true;

const unsigned int SHADOW_WIDTH = 2048, SHADOW_HEIGHT = 2048;
int WIDTH = 1000, HEIGHT = 1000;

float aspectRatio = 1.f;

glm::mat4 createPerspectiveMatrix()
{
	glm::mat4 perspectiveMatrix;
	float a1 = glm::min(aspectRatio, 1.f);
	float a2 = glm::min(1 / aspectRatio, 1.f);
	perspectiveMatrix = glm::mat4({
		1,0.,0.,0.,
		0.,aspectRatio,0.,0.,
		0.,0.,(far_plane + near_plane) / (near_plane - far_plane),2 * far_plane * near_plane / (near_plane - far_plane),
		0.,0.,-1.,0.,
		});

	perspectiveMatrix = glm::transpose(perspectiveMatrix);

	return perspectiveMatrix;
}


glm::mat4 createCameraMatrix()
{
	glm::vec3 cameraSide = glm::normalize(glm::cross(cameraDir, glm::vec3(0.f, 1.f, 0.f)));
	glm::vec3 cameraUp = glm::normalize(glm::cross(cameraSide, cameraDir));
	glm::mat4 cameraRotrationMatrix = glm::mat4({
		cameraSide.x,cameraSide.y,cameraSide.z,0,
		cameraUp.x,cameraUp.y,cameraUp.z ,0,
		-cameraDir.x,-cameraDir.y,-cameraDir.z,0,
		0.,0.,0.,1.,
		});
	cameraRotrationMatrix = glm::transpose(cameraRotrationMatrix);
	glm::mat4 cameraMatrix = cameraRotrationMatrix * glm::translate(-cameraPos);

	return cameraMatrix;
}



GLuint LoadCubemap(std::vector<std::string> faces) {
	GLuint textureID;
	glGenTextures(1, &textureID);
	glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);

	int width, height, nrChannels;
	for (GLuint i = 0; i < faces.size(); i++) {
		unsigned char* data = stbi_load(faces[i].c_str(), &width, &height, &nrChannels, 0);
		if (data) {
			glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
			stbi_image_free(data);
		}
		else {
			std::cerr << "Cubemap texture failed to load!" << std::endl;
			stbi_image_free(data);
			return 0;
		}
	}
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

	return textureID;
}

void drawSkybox(const glm::mat4& view, const glm::mat4& projection) {
	glDepthMask(GL_FALSE); // Wyłącz zapis do bufora głębokości
	glDisable(GL_DEPTH_TEST); // Wyłącz test głębokości

	glUseProgram(programSkybox);

	// Usuń translację z macierzy widoku
	glm::mat4 viewWithoutTranslation = glm::mat4(glm::mat3(view));

	// Przekaż macierz transformacji (połączona projekcja + widok)
	glm::mat4 transformation = projection * viewWithoutTranslation;

	// Przekaż macierz do shadera
	glUniformMatrix4fv(glGetUniformLocation(programSkybox, "transformation"), 1, GL_FALSE, &transformation[0][0]);

	// Aktywuj teksturę cubemapy
	glBindTexture(GL_TEXTURE_CUBE_MAP, skyboxTexture);

	// Narysuj skybox
	Core::DrawContext(skyboxContext);

	glDepthMask(GL_TRUE); // Włącz zapis do bufora głębokości
	glEnable(GL_DEPTH_TEST); // Włącz test głębokości

	glUseProgram(0); // Wyłącz program shaderów
}

namespace texture {
	GLuint earth;
	GLuint earthNormal;
	GLuint steel;
	GLuint steelNormal;
	GLuint dove;
	GLuint building;
	GLuint buildingNormal;
	GLuint bird;
}

float amountOfBoids = 50.f;
float neighborRadius = 1.0f;
float sightAngle = 120.0f;
float avoidBoids = 0.4f;
float avoidObstacles = 1.0f;
double lastTime = 0.0;
int attract = 0;
struct Boid {
	glm::vec3 position;   // Pozycja boida
	glm::vec3 velocity;   // Prędkość boida
	glm::vec3 acceleration; // Przyspieszenie boida
	glm::vec3 color;
	glm::vec3 minbox;
	glm::vec3 maxbox;
};
struct Obstacle {
	glm::vec3 position;  // Pozycja przeszkody
	float size;          // Rozmiar przeszkody (opcjonalnie)
	glm::vec3 minbox;
	glm::vec3 maxbox;
};
struct Building {
	glm::vec3 position;  // Pozycja przeszkody
	glm::vec3 size;          // Rozmiar przeszkody (opcjonalnie)
	glm::vec3 minbox;
	glm::vec3 maxbox;
};


std::vector<Boid> boids;
std::vector<Obstacle> obstacles;
std::vector<Building> buildings;

const int terrainWidth = 50;
const int terrainHeight = 50;
float terrain[terrainWidth][terrainHeight];
std::vector<float> vertices;
std::vector<unsigned int> indices;

//float frequencyValue = 0.1f;
 float frequencyValue = 0.5f; // -1.9

float Boundryfloat = 2.0f;
glm::vec3 minBoundary = glm::vec3(-10.0f, -Boundryfloat, -10.0f);
glm::vec3 maxBoundary = glm::vec3(10.0f, Boundryfloat, 10.0f);

void addBuilding(glm::vec3 buildPos, glm::vec3 buildSize = glm::vec3(0.5f, 1.0f, 1.0f)) {
	Building building;
	building.position = buildPos; // Pozycja przeszkody
	building.size = buildSize;          // Przykładowy rozmiar
	building.minbox = building.position - buildingboxsize * buildSize;
	building.maxbox = building.position + buildingboxsize * buildSize;
	buildings.push_back(building);
}


void initializeBoids(float numBoids, glm::vec3 color) {
	for (int i = 0; i < numBoids; i++) {
		Boid boid;
		boid.position = glm::vec3(rand() % 10 / 10.0f - 2.5f, rand() % 40 / 10.0f + 10.0f, rand() % 40 / 10.0f);
		boid.velocity = glm::vec3((rand() % 20 - 10) / 10.0f, (rand() % 20 - 10) / 10.0f, (rand() % 20 - 10) / 10.0f);
		boid.acceleration = glm::vec3(0.0f);
		boid.color = color;
		boids.push_back(boid);
	}
}

//obsluga wejscia
void processInput(GLFWwindow* window)
{
	if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, true);
	}

	float angleSpeed = 0.100f;
	float angleSpeed2 = 1.0f;
	float moveSpeed = 0.10f;
	if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
		cameraDir = glm::vec3(glm::eulerAngleY(angleSpeed) * glm::vec4(cameraDir, 0));
	}
	if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS) {
		cameraDir = glm::vec3(glm::eulerAngleY(-angleSpeed) * glm::vec4(cameraDir, 0));
	}
	if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS) {
		cameraPos += cameraDir * moveSpeed;
	}
	if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
		cameraPos -= cameraDir * moveSpeed;
	}
	if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
		cameraPos += glm::vec3(0, 1, 0) * moveSpeed;
	}
	if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
		cameraPos -= glm::vec3(0, 1, 0) * moveSpeed;
	}
	if (glfwGetKey(window, GLFW_KEY_O) == GLFW_PRESS) {
		Rotation1 -= angleSpeed2;
	}
	if (glfwGetKey(window, GLFW_KEY_P) == GLFW_PRESS) {
		Rotation2 -= angleSpeed2;
	}
}

void clearBoids() {
	boids.clear(); // Remove all boids from the vector
}

void RenderUI() {
	//std::cout << "RenderUI() executed" << std::endl;

	// Set position and size of the window (optional)
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(410, 200), ImGuiCond_Always);

	// Begin the ImGui window
	ImGui::Begin("UI Control Panel", NULL, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize);
	
	 // Static to keep it persistent across frames
	ImGui::SliderFloat("Boids:radius avoidBoid", &avoidBoids, 0.0f, 4.0f);
	ImGui::SliderFloat("Boids:radius avoidObstacle ", &avoidObstacles, 0.0f, 4.0f);
	ImGui::SliderFloat("Boids: Boundry ", &Boundryfloat, 0.0f, 10.0f);
	bool sliderChangedAmountBoids = ImGui::SliderFloat("Boids: number", &amountOfBoids, 0, 100);
	//bool sliderFrequencyChanged = = ImGui::SliderFloat("Terrain: Frequency", &amountOfBoids, 0, 100);

	// Check if the slider value has changed
	if (sliderChangedAmountBoids) {
		boids.clear();
		initializeBoids(amountOfBoids/2, glm::vec3(0.0, 1.0, 0.3));
		initializeBoids(amountOfBoids/2, glm::vec3(0.0, 0.0, 1.0));
	}

	ImGui::Dummy(ImVec2(5, 5));
	if (ImGui::Button(normalMappingEnabled ? "Disable Normal Mapping" : "Enable Normal Mapping")) {
		normalMappingEnabled = !normalMappingEnabled;
	}

	if (ImGui::Button("Spawn Obstacle")) {
		// Do something when the button is clicked
		Obstacle obstacle;
		obstacle.position = cameraPos + cameraDir; // Pozycja przeszkody
		obstacle.size = 0.1f;          // Przykładowy rozmiar
		obstacle.minbox = obstacle.position - obstacleboxsize;
		obstacle.maxbox = obstacle.position + obstacleboxsize;
		obstacles.push_back(obstacle);
	}

	if (ImGui::Button(isButtonClicked ? "Disable Attract" : "Enable Attract")) {
		attract = (attract + 1) % 2;
		isButtonClicked = !isButtonClicked;

		isPopupOpen = true;
		ImGui::OpenPopup(attract ? "Enable Popup Attract" : "Disable Popup Attract");
		popupOpenTime = std::chrono::steady_clock::now();
			
	}

	if (!isButtonClicked) {
		attract = 0;
	}

	if (isPopupOpen) {

		if (ImGui::BeginPopupModal("Enable Popup Attract", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
			ImGui::Text("Boid attraction is enabled!");

			auto now = std::chrono::steady_clock::now();
			if (std::chrono::duration<float>(now - popupOpenTime).count() >= 2.f) {
				isPopupOpen = false;
				ImGui::CloseCurrentPopup();
			}

			ImGui::EndPopup();  // End the popup
		}

		if (ImGui::BeginPopupModal("Disable Popup Attract", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
			ImGui::Text("Boid attraction is disabled!");

			auto now = std::chrono::steady_clock::now();
			if (std::chrono::duration<float>(now - popupOpenTime).count() >= 2.f) {
				isPopupOpen = false;
				ImGui::CloseCurrentPopup();
			}

			ImGui::EndPopup();  // End the popup
		}
	}

	ImGui::Dummy(ImVec2(5, 5));
	if (ImGui::Button(shadowMappingEnabled ? "Disable Shadow Mapping" : "Enable Shadow Mapping")) {
		shadowMappingEnabled = !shadowMappingEnabled;
	}
	// End the ImGui window5
	ImGui::End();
}



void InitImGui(GLFWwindow* window) {
	std::cout << "Initializing ImGui..." << std::endl;

	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	ImGuiIO& io = ImGui::GetIO();

	// Set the initial display size
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	io.DisplaySize = ImVec2((float)width, (float)height);

	io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
	std::cout << "Window size: " << io.DisplaySize.x << ", " << io.DisplaySize.y << std::endl;

	if (!ImGui_ImplGlfw_InitForOpenGL(window, true)) {
		std::cerr << "Failed to initialize ImGui GLFW backend!" << std::endl;
		return;
	}

	if (!ImGui_ImplOpenGL3_Init("#version 330")) {
		std::cerr << "Failed to initialize ImGui OpenGL backend!" << std::endl;
		return;
	}

	std::cout << "ImGui initialized successfully!" << std::endl;
}



// Function to generate terrain based on Perlin noise
void generateTerrain() {
	FastNoiseLite noise;
	noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
	noise.SetFrequency(frequencyValue); // Adjust zoom

	for (int x = 0; x < terrainWidth; x++) {
		for (int y = 0; y < terrainHeight; y++) {
			float value = noise.GetNoise((float)x * 0.05f, (float)y * 0.05f);
			terrain[x][y] = value - 1.9f;//2.1f; // Scale height for visibility
			//std::cout << "terrain[" << x << "][" << y << "] = " << terrain[x][y] << std::endl;
		}
	}
}

void createTerrainMesh(std::vector<float>& vertices, std::vector<unsigned int>& indices) {
	float halfWidth = terrainWidth / 2.0f;
	float halfHeight = terrainHeight / 2.0f;

	// Generate vertices
	for (int x = 0; x < terrainWidth; x++) {
		for (int y = 0; y < terrainHeight; y++) {
			vertices.push_back(x - halfWidth);   // X coordinate
			vertices.push_back(terrain[x][y] + 1.4f);  // Y (height)
			vertices.push_back(y - halfHeight); // Z coordinate

			float r = (terrain[x][y] + 1.4f) / 10.0f; // Normalize height for color (adjust as needed)
			float g = 0.5f;                // Fixed green value (adjust for different effects)
			float b = 1.0f - r;            // Inverse of red for variation

			vertices.push_back(r);
			vertices.push_back(g);
			vertices.push_back(b);

			// Texture Coordinates (Normalized)
			//float u = (float)x/halfWidth;  // U coordinate, normalized
			//float v = (float)y/halfHeight; // V coordinate, normalized
			//vertices.push_back(u);  // U
			//vertices.push_back(v);  // V
		}
	}

	// Generate indices for triangles
	for (int x = 0; x < terrainWidth - 1; x++) {
		for (int y = 0; y < terrainHeight - 1; y++) {
			int topLeft = x * terrainHeight + y;
			int topRight = (x + 1) * terrainHeight + y;
			int bottomLeft = x * terrainHeight + (y + 1);
			int bottomRight = (x + 1) * terrainHeight + (y + 1);

			// First Triangle (Top Left, Bottom Left, Top Right)
			indices.push_back(topLeft);
			indices.push_back(bottomLeft);
			indices.push_back(topRight);

			// Second Triangle (Top Right, Bottom Left, Bottom Right)
			indices.push_back(topRight);
			indices.push_back(bottomLeft);
			indices.push_back(bottomRight);
		}
	}
}

void setupBuildings() {
	int buildingSpacing = 6;  // Space between buildings
	int gridSize = 50;        // Grid size (same as terrain)

	for (int x = 0; x < gridSize; x += buildingSpacing) {
		for (int z = 0; z < gridSize; z += buildingSpacing) {
			float terrainHeight = terrain[x][z];  // Get terrain height at (x, z)

			// Ensure buildings are placed on relatively flat areas
			if (fabs(terrainHeight - terrain[x + 1][z]) < 0.5f &&
				fabs(terrainHeight - terrain[x][z + 1]) < 0.5f) {

				float height = 5.0f + (rand() % 10);  // Random building height (5-15 units)

				addBuilding(glm::vec3(x - gridSize / 2, terrainHeight + height / 2.0f, z - gridSize / 2), glm::vec3(2.0f, height, 2.0f));

			}
		}
	}
}

void initTerrain(const std::vector<float>& vertices, const std::vector<unsigned int>& indices) {
	glGenVertexArrays(1, &VAOTer);
	glGenBuffers(1, &VBOTer);
	glGenBuffers(1, &EBOTer);

	glBindVertexArray(VAOTer);

	glBindBuffer(GL_ARRAY_BUFFER, VBOTer);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBOTer);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	// Position attribute
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	// Color attribute
	glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
	glEnableVertexAttribArray(1);

	glBindVertexArray(0);
	

}

// Function to render the terrain
void renderTerrain() {
	//glUseProgram(programTerrain);

	//glm::mat4 model = glm::mat4(1.0f);
	//glm::mat4 view = createCameraMatrix();
	//glm::mat4 projection = createPerspectiveMatrix();

	//glUniformMatrix4fv(glGetUniformLocation(programTerrain, "model"), 1, GL_FALSE, glm::value_ptr(model));
	//glUniformMatrix4fv(glGetUniformLocation(programTerrain, "view"), 1, GL_FALSE, glm::value_ptr(view));
	//glUniformMatrix4fv(glGetUniformLocation(programTerrain, "projection"), 1, GL_FALSE, glm::value_ptr(projection));

	glBindVertexArray(VAOTer);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}

bool checkCollision(const glm::vec3& min1, const glm::vec3& max1, const glm::vec3& min2, const glm::vec3& max2) {
	// Sprawdzamy, czy zakresy na każdej z osi się nie nakładają
	if (max1.x < min2.x || min1.x > max2.x) return false; // Oś X
	if (max1.y < min2.y || min1.y > max2.y) return false; // Oś Y
	if (max1.z < min2.z || min1.z > max2.z) return false; // Oś Z

	// Jeśli żadna z osi nie jest rozłączna, to kolizja występuje
	return true;
}

void checkCollisionBoid(Boid& boid) {
	for (const auto& obstacle : obstacles) {

		if (checkCollision(boid.minbox, boid.maxbox, obstacle.minbox, obstacle.maxbox)) {
			boid.color = glm::vec3(1.0, 0.0, 0.0);
		}
	}
	for (const auto& building : buildings) {

		if (checkCollision(boid.minbox, boid.maxbox, building.minbox, building.maxbox)) {
			boid.color = glm::vec3(1.0, 0.0, 0.0);
		}
	}
}

bool insight(const Boid& main, glm::vec3 position, float radius) {
	glm::vec3 directionToBoid2 = position - main.position;

	float dotProduct = glm::dot(main.velocity, directionToBoid2);

	// Obliczamy długości wektorów
	float magnitudeVelocity = glm::length(main.velocity);
	float magnitudeDirection = glm::length(directionToBoid2);

	// Obliczamy cos kąta
	float cosTheta = dotProduct / (magnitudeVelocity * magnitudeDirection);

	// Zabezpieczenie przed przekroczeniem zakresu, aby uniknąć błędów numerycznych
	cosTheta = glm::clamp(cosTheta, -1.0f, 1.0f);

	// Obliczamy kąt w radianach
	float angle = acos(cosTheta);

	if (glm::degrees(angle) > sightAngle) {
		return false; // Boid nie widzi drugiego boida, ponieważ jest poza polem widzenia
	}
	float distance = glm::distance(main.position, position);

	if (distance > 0 && distance < radius) {
		return true;
	}

	return false;
}

void checkPosition(Boid& boid, glm::vec3 minBoundary, glm::vec3 maxBoundary) {
	// Sprawdź granice w osi X
	if (boid.position.x > maxBoundary.x) {
		boid.position.x = maxBoundary.x;  // Ustaw pozycję na granicy
		boid.velocity.x = -boid.velocity.x;  // Odbij prędkość
	}
	if (boid.position.x < minBoundary.x) {
		boid.position.x = minBoundary.x;
		boid.velocity.x = -boid.velocity.x;
	}

	// Sprawdź granice w osi Z
	if (boid.position.z > maxBoundary.z) {
		boid.position.z = maxBoundary.z;
		boid.velocity.z = -boid.velocity.z;
	}
	if (boid.position.z < minBoundary.z) {
		boid.position.z = minBoundary.z;
		boid.velocity.z = -boid.velocity.z;
	}

	// Sprawdź granice w osi Y
	if (boid.position.y > maxBoundary.y) {
		boid.position.y = maxBoundary.y;  // Ustaw pozycję na granicy
		boid.velocity.y = -boid.velocity.y;  // Odbij prędkość
	}
	if (boid.position.y < minBoundary.y) {
		boid.position.y = minBoundary.y;
		boid.velocity.y = -boid.velocity.y;
	}
}



glm::vec3 cohesion(Boid& boid, const std::vector<Boid>& boids, float neighborRadius) {
	glm::vec3 centerOfMass(0.0f);
	int count = 0;

	for (const auto& other : boids) {
		if (insight(boid, other.position, neighborRadius) && boid.color == other.color) {
			centerOfMass += other.position;
			count++;
		}
	}

	if (count > 0) {
		centerOfMass /= count;
		return glm::normalize(centerOfMass - boid.position);
	}
	return glm::vec3(0.0f);
}

glm::vec3 separation(Boid& boid, const std::vector<Boid>& boids, float avoidBoids) {
	glm::vec3 avoid(0.0f);  // Siła unikania
	int count = 0;

	// Przeszukaj wszystkie boidy
	for (const auto& other : boids) {
		// Upewnijmy się, że boidy są w zasięgu widzenia i w obrębie unikania
		if (insight(boid, other.position, avoidBoids) && boid.color == other.color) {
			float distance = glm::distance(boid.position, other.position);
			if (distance < avoidBoids) {
				glm::vec3 direction = boid.position - other.position;
				avoid += glm::normalize(direction) / (distance);  // Normalizuj wektor, aby uwzględnić odległość
				count++;
			}
		}
	}


	if (count > 0) {
		avoid /= count;  // Przeciętna siła unikania
	}

	// Jeśli siła unikania jest za mała, zwróć zero, aby uniknąć błędów
	if (glm::length(avoid) > 0.0f) {
		return glm::normalize(avoid);  // Znormalizuj, aby siła była stała
	}

	return glm::vec3(0.0f);  // W przeciwnym razie, nie generuj żadnej siły
}

glm::vec3 separationObstacles(Boid& boid, const std::vector<Obstacle>& obstacles, float avoidObstacles) {
	glm::vec3 avoid(0.0f);  // Siła unikania
	int count = 0;


	for (const auto& other : obstacles) {
		// Upewnijmy się, że boidy są w zasięgu widzenia i w obrębie unikania
		if (insight(boid, other.position, avoidObstacles)) {
			float distance = glm::distance(boid.position, other.position);
			if (distance < avoidBoids) {
				glm::vec3 direction = boid.position - other.position;
				avoid += glm::normalize(direction) / (distance);  // Normalizuj wektor, aby uwzględnić odległość
				count++;
			}
		}
	}


	if (count > 0) {
		avoid /= count;  // Przeciętna siła unikania
	}

	// Jeśli siła unikania jest za mała, zwróć zero, aby uniknąć błędów
	if (glm::length(avoid) > 0.0f) {
		return glm::normalize(avoid);  // Znormalizuj, aby siła była stała
	}

	return glm::vec3(0.0f);
}

glm::vec3 separationBuildings(Boid& boid, const std::vector<Building>& buildings) {
	glm::vec3 avoid(0.0f);  // Siła unikania
	int count = 0;

	for (const auto& building : buildings) {
		if (checkCollision(boid.minbox, boid.maxbox, building.minbox, building.maxbox)) {
			//std::cout << "Boid Position: ("
				//<< boid.position.x << ", "
				//<< boid.position.y << ", "
				//<< boid.position.z << ")" << std::endl;

			// Jeśli kolizja zachodzi, dodaj wektor unikania
			glm::vec3 closestPoint = glm::clamp(boid.position, building.minbox, building.maxbox);
			glm::vec3 direction = boid.position - closestPoint;
			float distance = glm::length(direction);

			// Upewnij się, że wektor ma długość > 0
			if (distance > 0.0f) {
				avoid += glm::normalize(direction) / distance;
				count++;
			}
		}
	}

	if (count > 0) {
		avoid /= count;  // Średnia siła unikania
	}

	// Normalizacja siły unikania, jeśli istnieje
	if (glm::length(avoid) > 0.0f) {
		return glm::normalize(avoid);
	}

	return glm::vec3(0.0f);  // Brak siły unikania
}

glm::vec3 alignment(Boid& boid, const std::vector<Boid>& boids, float neighborRadius) {
	glm::vec3 averageVelocity(0.0f);
	int count = 0;

	for (const auto& other : boids) {
		if (insight(boid, other.position, neighborRadius) && boid.color == other.color) {
			averageVelocity += other.velocity;
			count++;
		}
	}

	if (count > 0) {
		averageVelocity /= count;
		return glm::normalize(averageVelocity);
	}
	return glm::vec3(0.0f);
}

glm::vec3 attraction(Boid& boid) {
	glm::vec3 force = glm::vec3(0.0f);
	if (attract == 1) {
		force = cameraPos + cameraDir - boid.position;

		// Oblicz odległość między boidem a statkiem
		float distance = glm::length(force);

		// Jeśli odległość jest bardzo mała, unikamy dzielenia przez zero
		if (distance < 0.01f) {
			return glm::vec3(0.0f); // Siła zerowa, bo są bardzo blisko siebie
		}

		// Normalizuj wektor siły i dziel przez odległość (siła maleje z odległością)
		force = glm::normalize(force) / distance;
	}
	return force;
}

glm::vec3 limitSpeed(const glm::vec3& vector, float maxLength) {
	float length = glm::length(vector);
	if (length > maxLength) {
		return glm::normalize(vector) * maxLength;
	}
	return vector;
}

void updateboxBoid(Boid& boid, glm::vec3 boxsize) {
	boid.minbox = boid.position - boxsize;
	boid.maxbox = boid.position + boxsize;
}


void updateBoids(float deltaTime, float neighborRadius, float avoidBoids) {

	for (auto& boid : boids) {
		glm::vec3 forceCohesion = cohesion(boid, boids, neighborRadius);
		glm::vec3 forceSeparation = separation(boid, boids, avoidBoids);
		glm::vec3 forceSeparationObstacles = separationObstacles(boid, obstacles, avoidObstacles);
		glm::vec3 forceSeparationBuildings = separationBuildings(boid, buildings);
		glm::vec3 forceAlignment = alignment(boid, boids, neighborRadius);
		glm::vec3 forceattraction = attraction(boid);

		boid.acceleration = forceAlignment + forceCohesion + forceSeparation + forceSeparationObstacles * 5 + forceSeparationBuildings * 10 + forceattraction;
		boid.velocity += boid.acceleration * deltaTime;
		boid.velocity = limitSpeed(boid.velocity, maxSpeed);
		boid.position += boid.velocity * deltaTime;

		
		checkPosition(boid, minBoundary, maxBoundary);

		updateboxBoid(boid, boxsize);
		checkCollisionBoid(boid);
	}
}



void drawObjectTexture(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint texture) {

	glUseProgram(programTex);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(programTex, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(programTex, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(programTex, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	glUniform3f(glGetUniformLocation(programTex, "lightDir"), lightDir.x, lightDir.y, lightDir.z);
	Core::SetActiveTexture(texture, "colorTexture", programTex, 0);

	Core::DrawContext(context);
}

void drawObjectTextureShadow(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint texture) {

	glUseProgram(programTexShadow);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(programTexShadow, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(programTexShadow, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(programTexShadow, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	glUniform3f(glGetUniformLocation(programTexShadow, "lightDir"), lightDir.x, lightDir.y, lightDir.z);
	Core::SetActiveTexture(texture, "colorTexture", programTexShadow, 0);

	glUniformMatrix4fv(glGetUniformLocation(programTexShadow, "LightVP"), 1, GL_FALSE, (float*)&lightVP);
	Core::SetActiveTexture(depthMap, "depthMap", programTexShadow, 1);

	Core::DrawContext(context);
}

void drawObjectNormal(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint textureID, GLuint normalmapId) {
	glUseProgram(programNormal);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(programNormal, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(programNormal, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(programNormal, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	glUniform3f(glGetUniformLocation(programNormal, "lightDir"), lightDir.x, lightDir.y, lightDir.z);

	Core::SetActiveTexture(textureID, "colorTexture", programNormal, 0);
	Core::SetActiveTexture(normalmapId, "normalSampler", programNormal, 1);

	Core::DrawContext(context);
}

void drawObjectNormalShadow(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint textureID, GLuint normalmapId) {
	glUseProgram(programNormalShadow);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(programNormalShadow, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(programNormalShadow, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(programNormalShadow, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	glUniform3f(glGetUniformLocation(programNormalShadow, "lightDir"), lightDir.x, lightDir.y, lightDir.z);

	Core::SetActiveTexture(textureID, "colorTexture", programNormalShadow, 0);
	Core::SetActiveTexture(normalmapId, "normalSampler", programNormalShadow, 1);

	glUniformMatrix4fv(glGetUniformLocation(programNormalShadow, "LightVP"), 1, GL_FALSE, (float*)&lightVP);
	Core::SetActiveTexture(depthMap, "depthMap", programNormalShadow, 2);

	Core::DrawContext(context);
}

void drawObjectDepth(Core::RenderContext context, glm::mat4 viewProjectionMatrix, glm::mat4 modelMatrix) {

	glUniformMatrix4fv(glGetUniformLocation(shadowShaderProgram, "viewProjectionMatrix"), 1, GL_FALSE, (float*)&viewProjectionMatrix);
	glUniformMatrix4fv(glGetUniformLocation(shadowShaderProgram, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	Core::DrawContext(context);
}

void drawBoids(bool shadow, glm::mat4 lightViewProjectionMatrix) {
	for (const auto& boid : boids) {
		glm::vec3 forward = glm::normalize(boid.velocity);
		glm::vec3 Up = glm::vec3(0, 1, 0);
		glm::vec3 right = glm::normalize(glm::cross(Up, forward));

		glm::vec3 boidUp = glm::cross(forward, right);

		glm::mat4 rotationMatrix = glm::mat4({
		right.x,right.y,right.z,0.,
		boidUp.x,boidUp.y,boidUp.z,0.,
		forward.x,forward.y,forward.z,0.,
		0.,0.,0.,1.,
			});
		glm::mat4 translationMatrix = glm::translate(glm::mat4(1.0f), boid.position);


		glm::mat4 modelMatrix = translationMatrix * rotationMatrix;

		if (shadow) {
			drawObjectDepth(coneContext, lightViewProjectionMatrix, modelMatrix);
		}
		else {
			if (shadowMappingEnabled) {
				drawObjectTextureShadow(coneContext, modelMatrix, texture::bird);
			}
			else {
				drawObjectTexture(coneContext, modelMatrix, texture::bird);
			}
		}
	}
}

void drawObstacles(bool shadow, glm::mat4 lightViewProjectionMatrix) {
	for (const auto& obstacle : obstacles) {
		// Najpierw translacja
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), obstacle.position);

		// Następnie skalowanie
		modelMatrix = glm::scale(modelMatrix, glm::vec3(obstacle.size, obstacle.size, obstacle.size));

		if (shadow) {
			drawObjectDepth(sphereContext, lightViewProjectionMatrix, modelMatrix);
		}
		else {
			if (normalMappingEnabled) {
				if (shadowMappingEnabled) {
					drawObjectNormalShadow(sphereContext, modelMatrix, texture::earth, texture::earthNormal);
				}
				else {
					drawObjectNormal(sphereContext, modelMatrix, texture::earth, texture::earthNormal);
				}
			}
			else {
				if (shadowMappingEnabled) {
					drawObjectTextureShadow(sphereContext, modelMatrix, texture::earth);
				}
				else {
					drawObjectTexture(sphereContext, modelMatrix, texture::earth);
				}
			}
		}
	}
}

void drawBuildings(bool shadow, glm::mat4 lightViewProjectionMatrix) {
	for (const auto& building : buildings) {
		// Najpierw translacja
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), building.position);

		// Następnie skalowanie
		modelMatrix = glm::scale(modelMatrix, glm::vec3(building.size.x, building.size.y, building.size.z));

		if (shadow) {
			drawObjectDepth(buildingContext, lightViewProjectionMatrix, modelMatrix);
		}
		else {
			if (normalMappingEnabled) {
				if (shadowMappingEnabled) {
					drawObjectNormalShadow(buildingContext, modelMatrix, texture::building, texture::buildingNormal);
				}
				else {
					drawObjectNormal(buildingContext, modelMatrix, texture::building, texture::buildingNormal);
				}
			}
			else {
				if (shadowMappingEnabled) {
					drawObjectTextureShadow(buildingContext, modelMatrix, texture::building);
				}
				else {
					drawObjectTexture(buildingContext, modelMatrix, texture::building);
				}
			}
		}
	}
}

void makescene(bool shadow, glm::mat4 lightViewProjectionMatrix) {

	drawBuildings(shadow, lightViewProjectionMatrix);
}

void renderShadow() {
	float time = glfwGetTime();
	glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);

	glUseProgram(shadowShaderProgram);
	glViewport(0, 0, SHADOW_WIDTH, SHADOW_HEIGHT);
	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO); 
	glClear(GL_DEPTH_BUFFER_BIT);

	glUseProgram(shadowShaderProgram);

	drawBoids(true, lightVP);
	makescene(true, lightVP);
	drawObstacles(true, lightVP);

	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glViewport(0, 0, WIDTH, HEIGHT);
}

void initDepthMap()
{
	glGenFramebuffers(1, &depthMapFBO);

	glGenTextures(1, &depthMap);
	glBindTexture(GL_TEXTURE_2D, depthMap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT,
		SHADOW_WIDTH, SHADOW_HEIGHT, 0, GL_DEPTH_COMPONENT, GL_FLOAT, NULL);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

	float borderColor[] = { 1.0f, 1.0f, 1.0f, 1.0f };
	glTexParameterfv(GL_TEXTURE_2D, GL_TEXTURE_BORDER_COLOR, borderColor);

	glBindFramebuffer(GL_FRAMEBUFFER, depthMapFBO);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depthMap, 0);
	glDrawBuffer(GL_NONE);
	glReadBuffer(GL_NONE);
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void renderScene(GLFWwindow* window)
{
	glClearColor(0.0f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 view = createCameraMatrix();
	glm::mat4 projection = createPerspectiveMatrix();
	drawSkybox(view, projection);

	glUseProgram(program);

	if (shadowMappingEnabled) {
		glCullFace(GL_FRONT);
		renderShadow();
		glCullFace(GL_BACK);
	}

	double currentTime = glfwGetTime();
	double time = currentTime - lastTime;
	updateBoids(time, neighborRadius, avoidBoids);
	lastTime = currentTime;

	drawBoids(false, glm::mat4(0.f));
	makescene(false, glm::mat4(0.f));
	renderTerrain();
	drawObstacles(false, glm::mat4(0.f));
	
	glm::mat4 cameraMatrix = createCameraMatrix();

	ImGui_ImplOpenGL3_NewFrame();
	ImGui_ImplGlfw_NewFrame();
	ImGui::NewFrame();

	// Render UI
	RenderUI();

	ImGui::Render();
	ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

	glUseProgram(0);

	glfwSwapBuffers(window);
}

void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	aspectRatio = width / float(height);
	glViewport(0, 0, width, height);
	WIDTH = width;
	HEIGHT = height;
}

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 texCoords;
};

void loadModelToContext(std::string path, Core::RenderContext& context)
{
	Assimp::Importer import;
	const aiScene* scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_CalcTangentSpace);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
	{
		std::cout << "ERROR::ASSIMP::" << import.GetErrorString() << std::endl;
		return;
	}
	context.initFromAssimpMesh(scene->mMeshes[0]);
}

void init(GLFWwindow* window)
{
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	initDepthMap();
	glEnable(GL_DEPTH_TEST);
	program = shaderLoader.CreateProgram("shaders/shader.vert", "shaders/shader.frag");
	programTex = shaderLoader.CreateProgram("shaders/shader_tex.vert", "shaders/shader_tex.frag");
	shadowShaderProgram = shaderLoader.CreateProgram("shaders/shader_shadow.vert", "shaders/shader_shadow.frag");
	programNormal = shaderLoader.CreateProgram("shaders/shader_normal.vert", "shaders/shader_normal.frag");
	programNormalShadow = shaderLoader.CreateProgram("shaders/shader_normal_shadow.vert", "shaders/shader_normal_shadow.frag");
	programTexShadow = shaderLoader.CreateProgram("shaders/shader_tex_shadow.vert", "shaders/shader_tex_shadow.frag");
	programTerrain = shaderLoader.CreateProgram("shaders/shader_terrain.vert", "shaders/shader_terrain.frag");
	loadModelToContext("./models/dove.obj", coneContext);
	loadModelToContext("./models/sphere.obj", sphereContext);
	loadModelToContext("./models/cuboid.obj", buildingContext);
	loadModelToContext("./models/spaceship.obj", shipContext);
	initializeBoids(amountOfBoids, glm::vec3(0.0, 1.0, 0.3));
	initializeBoids(amountOfBoids, glm::vec3(0.0, 0.0, 1.0));

	texture::earth = Core::LoadTexture("./textures/earth.png");
	texture::earthNormal = Core::LoadTexture("./textures/earth_normalmap.png");

	texture::steel = Core::LoadTexture("./textures/Steel_S.jpg");
	texture::steelNormal = Core::LoadTexture("./textures/Steel_N.jpg");

	texture::dove = Core::LoadTexture("./textures/DOVE.jpg");

	texture::building = Core::LoadTexture("./textures/building.jpg");
	texture::buildingNormal = Core::LoadTexture("./textures/buildingNormal.jpg");
	texture::bird = Core::LoadTexture("./textures/ptak.jpg");
	loadModelToContext("./models/cube.obj", skyboxContext);

	std::vector<std::string> skyboxFaces = {
		"textures/skybox/right_c.jpg",
		"textures/skybox/left.jpg",
		"textures/skybox/top.jpg",
		"textures/skybox/bottom.jpg",
		"textures/skybox/back.jpg",
		"textures/skybox/front.jpg"
	};

	skyboxTexture = LoadCubemap(skyboxFaces);
	programSkybox = shaderLoader.CreateProgram("shaders/shader_skybox.vert", "shaders/shader_skybox.frag");
}

void shutdown(GLFWwindow* window)
{
	shaderLoader.DeleteProgram(program);
}



void renderLoop(GLFWwindow* window) {
	InitImGui(window);

	generateTerrain();
	setupBuildings();
	createTerrainMesh(vertices, indices);
	initTerrain(vertices, indices);

	while (!glfwWindowShouldClose(window))
	{
		processInput(window);

		// Clear OpenGL buffers
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

		// Start ImGui frame

		GLuint cameraPosLocation = glGetUniformLocation(programNormal, "cameraPos");
		glUseProgram(programNormal);
		glUniform3fv(cameraPosLocation, 1, glm::value_ptr(cameraPos));
		glUseProgram(0);

		renderScene(window);

		glfwPollEvents();
	}
}