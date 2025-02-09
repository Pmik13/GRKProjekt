#include "glew.h"
#include <GLFW/glfw3.h>
#include "glm.hpp"
#include "ext.hpp"
#include <iostream>
#include <cmath>

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
#include <unordered_map>

GLuint program;
GLuint programTex;
GLuint programNormal;
Core::Shader_Loader shaderLoader;
GLuint shadowShaderProgram;
GLuint VAO, VBO;
GLuint depthMapFBO;
GLuint depthMap;
GLuint skyboxTexture;
GLuint programSkybox;
Core::RenderContext skyboxContext;
Core::RenderContext buildingContext;
Core::RenderContext coneContext;
Core::RenderContext sphereContext;
float Rotation1 = 0.0f;
float Rotation2 = 0.0f;
float Rotation3 = 0.0f;
float cameraAngle = 0;
float maxSpeed = 1.0f;

glm::vec3 obstacleboxsize = glm::vec3(0.05f, 0.05f, 0.05f);
glm::vec3 buildingboxsize = glm::vec3(0.5f, 0.5f, 0.5f);
glm::vec3 boxsize = glm::vec3(0.05f, 0.05f, 0.05f);
glm::vec3 cameraPos = glm::vec3(-4.f, 0.0f, 0.0f);
glm::vec3 cameraDir = glm::vec3(1.f, 0.f, 0.f);;
glm::vec3 spaceshipPos = glm::vec3(-4.f, 0, 0);
glm::vec3 spaceshipDir = glm::vec3(1.f, 0.f, 0.f);

float near_plane = 0.05;
float far_plane = 20.;
glm::vec3 lightPos = glm::vec3(0.f, 1.f, 0.f);

bool shadowMappingEnabled = false;
bool normalMappingEnabled = true;

float aspectRatio = 1.f;

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
}

float amountOfBoids = 1.f;
float neighborRadius = 1.0f;
float sightAngle = 120.0f;
float avoidBoids = 0.4f;
float avoidObstacles = 1.0f;
double lastTime = 0.0;
int attract = 0;

std::vector<glm::vec3> compute12DOP(const std::vector<glm::vec3>& vertices, const std::vector<unsigned int>& indices) {
	std::vector<glm::vec3> DOP(12);  // Wyniki
	std::vector<glm::vec3> directions = {
		glm::vec3(1, 0, 0), glm::vec3(-1, 0, 0),
		glm::vec3(0, 1, 0), glm::vec3(0, -1, 0),
		glm::vec3(0, 0, 1), glm::vec3(0, 0, -1),
		glm::normalize(glm::vec3(1, 1, 0)), glm::normalize(glm::vec3(-1, 1, 0)),
		glm::normalize(glm::vec3(1, -1, 0)), glm::normalize(glm::vec3(-1, -1, 0)),
		glm::normalize(glm::vec3(1, 0, 1)), glm::normalize(glm::vec3(-1, 0, 1))
	};

	// Iteruj przez wszystkie kierunki
	for (int i = 0; i < 12; i++) {
		float minProj = glm::dot(vertices[indices[0]], directions[i]);  // Projektuj pierwszy wierzchołek
		float maxProj = minProj;  // Inicjalizuj max i min dla kierunku

		// Iteruj przez wszystkie indeksy wierzchołków
		for (size_t j = 0; j < indices.size(); j++) {
			unsigned int idx = indices[j];
			glm::vec3 vertex = vertices[idx];  // Pobierz wierzchołek za pomocą indeksu
			float proj = glm::dot(vertex, directions[i]);

			// Zaktualizuj max i min
			if (proj > maxProj) maxProj = proj;
			if (proj < minProj) minProj = proj;
		}

		// Przypisz wynik dla danego kierunku
		DOP[i] = directions[i] * (maxProj - minProj);
	}

	return DOP;
}


struct Boid {
	glm::vec3 position;   // Pozycja boida
	glm::vec3 velocity;   // Prędkość boida
	glm::vec3 acceleration; // Przyspieszenie boida
	glm::vec3 color;
	std::vector<glm::vec3> vertices;
	std::vector<unsigned int> indices;
	std::vector<glm::vec3> DOP;

	void setKDOP() {
		DOP = compute12DOP(vertices, indices);
	}
};

struct Obstacle {
	glm::vec3 position;  // Pozycja przeszkody
	float size;          // Rozmiar przeszkody (opcjonalnie)
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> DOP;
	std::vector<unsigned int> indices;

	void setKDOP() {
		DOP = compute12DOP(vertices, indices);
	}
};

struct Building {
	glm::vec3 position;
	glm::vec3 size;
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec3> DOP;
	std::vector<unsigned int> indices;

	void setKDOP() {
		DOP = compute12DOP(vertices, indices);
	}
};

std::vector<Boid> boids;
std::vector<Obstacle> obstacles;
std::vector<Building> buildings;

const int terrainWidth = 50;
const int terrainHeight = 50;
float terrain[terrainWidth][terrainHeight];
std::vector<float> vertices;
std::vector<unsigned int> indices;

float frequencyValue = 0.9f;

float Boundryfloat = 2.0f;
glm::vec3 minBoundary = glm::vec3(-Boundryfloat, -Boundryfloat, -Boundryfloat);
glm::vec3 maxBoundary = glm::vec3(Boundryfloat, Boundryfloat, Boundryfloat);

void addBuilding(glm::vec3 buildPos, glm::vec3 buildSize) {
	Building building;
	building.position = buildPos;
	building.size = buildSize;
	building.vertices = buildingContext.getVertices();
	building.indices = buildingContext.getIndices();
	building.setKDOP();

	buildings.push_back(building);
}

void initializeBoids(float numBoids, glm::vec3 color) {
	for (int i = 0; i < numBoids; i++) {
		Boid boid;
		boid.position = glm::vec3(rand() % 5 / 10.0f - 2.5f, rand() % 10 / 10.0f, rand() % 10 / 10.0f);
		boid.velocity = glm::vec3((rand() % 20 - 10) / 10.0f, (rand() % 20 - 10) / 10.0f, (rand() % 20 - 10) / 10.0f);
		boid.acceleration = glm::vec3(0.0f);
		boid.color = color;
		boid.indices = coneContext.getIndices();
		boid.vertices = coneContext.getVertices();

		boid.setKDOP();
		boids.push_back(boid);
	}
}

void clearBoids() {
	boids.clear(); // Remove all boids from the vector
}

void RenderUI() {
	//std::cout << "RenderUI() executed" << std::endl;

	// Set position and size of the window (optional)
	ImGui::SetNextWindowPos(ImVec2(10, 10), ImGuiCond_Always);
	ImGui::SetNextWindowSize(ImVec2(410, 160), ImGuiCond_Always);

	// Begin the ImGui window
	ImGui::Begin("UI Control Panel", NULL, ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize);
	
	 // Static to keep it persistent across frames
	ImGui::SliderFloat("Boids:radius avoidBoid", &avoidBoids, 0.0f, 4.0f);
	ImGui::SliderFloat("Boids:radius avoidObstacle ", &avoidObstacles, 0.0f, 4.0f);
	ImGui::SliderFloat("Boids: Boundry ", &Boundryfloat, 0.0f, 10.0f);
	bool sliderChanged = ImGui::SliderFloat("Boids: number", &amountOfBoids, 0, 100);


	// Check if the slider value has changed
	if (sliderChanged) {
		boids.clear();
		initializeBoids(amountOfBoids/2, glm::vec3(0.0, 1.0, 0.3));
		initializeBoids(amountOfBoids/2, glm::vec3(0.0, 0.0, 1.0));
	}

	ImGui::Dummy(ImVec2(5, 5));
	if (ImGui::Button(normalMappingEnabled ? "Disable Normal Mapping" : "Enable Normal Mapping")) {
		normalMappingEnabled = !normalMappingEnabled;
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
			float value = noise.GetNoise((float)x * 0.1f, (float)y * 0.1f);
			terrain[x][y] = value - 2.0f; // Scale height for visibility
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
	int buildingSpacing = 5;  // Space between buildings
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

// Function to render the terrain
void renderTerrain(const std::vector<float>& vertices, const std::vector<unsigned int>& indices) {
	// Create and bind VAO, VBO, and EBO
	unsigned int VAO, VBO, EBO;
	glGenVertexArrays(1, &VAO);
	glGenBuffers(1, &VBO);
	glGenBuffers(1, &EBO);

	glBindVertexArray(VAO);

	glBindBuffer(GL_ARRAY_BUFFER, VBO);
	glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

	// Define vertex attribute pointers
	glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
	glEnableVertexAttribArray(0);

	glBindBuffer(GL_ARRAY_BUFFER, 0);
	glBindVertexArray(0);

	// Use your shader program
	// glUseProgram(shaderProgram);

	// Bind the VAO and draw the terrain
	glBindVertexArray(VAO);
	glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
	glBindVertexArray(0);
}

bool checkKDOPCollision(const std::vector<glm::vec3>& DOP1, const std::vector<glm::vec3>& DOP2) {
	// Iterujemy przez 12 osi (w przypadku 12DOP)
	for (int i = 0; i < 12; i++) {
		// Rzutowanie DOP1
		float minProj1 = glm::dot(DOP1[0], DOP1[i]);
		float maxProj1 = minProj1;

		// Zamiast iterować przez wszystkie wierzchołki DOP1, tylko wektory, które tworzą 12DOP
		for (int j = 1; j < DOP1.size(); j++) {
			float proj = glm::dot(DOP1[j], DOP1[i]);
			if (proj > maxProj1) maxProj1 = proj;
			if (proj < minProj1) minProj1 = proj;
		}

		// Rzutowanie DOP2
		float minProj2 = glm::dot(DOP2[0], DOP2[i]);
		float maxProj2 = minProj2;

		// Zamiast iterować przez wszystkie wierzchołki DOP2, tylko wektory, które tworzą 12DOP
		for (int j = 1; j < DOP2.size(); j++) {
			float proj = glm::dot(DOP2[j], DOP2[i]);
			if (proj > maxProj2) maxProj2 = proj;
			if (proj < minProj2) minProj2 = proj;
		}

		// Jeśli przedziały rzutów się nie zachodzą, to nie ma kolizji
		if (maxProj1 < minProj2 || maxProj2 < minProj1) {
			return false;
		}
	}
	return true;
}

void updateKDOPBoid(Boid& boid) {
	glm::vec3 delta = boid.velocity;

	for (int i = 0; i < 12; i++) {
		float deltaProj = glm::dot(delta, glm::normalize(boid.DOP[i]));
		boid.DOP[i] += glm::normalize(boid.DOP[i]) * deltaProj;
	}
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

	for (const auto& obstacle : obstacles) {
		if (checkKDOPCollision(boid.DOP, obstacle.DOP)) { // Sprawdzenie KDOP zamiast pozycji
			glm::vec3 direction = boid.position - obstacle.position;
			float distance = glm::length(direction);

			if (distance > 0.0f && distance < avoidObstacles) {
				avoid += glm::normalize(direction) / (distance + 0.1f); // Mniejszy wpływ dalekich obiektów
				count++;
			}
		}
	}

	if (count > 0) {
		avoid /= count;
		avoid *= 0.8f; // Zmniejszenie wpływu unikania, żeby boidy nie "drżały"
	}

	if (glm::length(avoid) > 0.0f) {
		return glm::normalize(avoid);  // Znormalizuj, aby siła była stała
	}

	return avoid;
}

glm::vec3 separationBuildings(Boid& boid, const std::vector<Building>& buildings) {
	glm::vec3 avoid(0.0f);
	int count = 0;

	for (const auto& building : buildings) {
		if (checkKDOPCollision(boid.DOP, building.DOP)) {
			glm::vec3 closestPoint = glm::clamp(boid.position,
				building.position - building.size / 2.0f,
				building.position + building.size / 2.0f);
			glm::vec3 direction = boid.position - closestPoint;
			float distance = glm::length(direction);

			if (distance > 0.0f) {
				avoid += glm::normalize(direction) / (distance + 0.2f);
				count++;
			}
		}
	}

	return (count > 0) ? glm::normalize(avoid) * 0.4f : glm::vec3(0.0f);
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

std::unordered_map<int, std::vector<Boid*>> spatialGrid;
float gridSize = 1.0f; // Rozmiar komórki siatki

int getGridKey(glm::vec3 position) {
	return (int)(position.x / gridSize) * 73856093
		^ (int)(position.y / gridSize) * 19349663
		^ (int)(position.z / gridSize) * 83492791;
}

void updateBoids(float deltaTime, float neighborRadius, float avoidBoids) {
	spatialGrid.clear(); // Reset siatki przestrzennej

	for (auto& boid : boids) {
		int gridKey = getGridKey(boid.position);
		spatialGrid[gridKey].push_back(&boid);
	}

	for (auto& boid : boids) {
		glm::vec3 forceCohesion(0.0f), forceSeparation(0.0f), forceAlignment(0.0f);
		int count = 0;

		// Sprawdź aktualną komórkę + sąsiadujące
		for (int dx = -1; dx <= 1; dx++) {
			for (int dy = -1; dy <= 1; dy++) {
				for (int dz = -1; dz <= 1; dz++) {
					int neighborKey = getGridKey(boid.position + glm::vec3(dx * gridSize, dy * gridSize, dz * gridSize));
					if (spatialGrid.find(neighborKey) != spatialGrid.end()) {
						for (auto& other : spatialGrid[neighborKey]) {
							if (other == &boid) continue;
							forceCohesion += cohesion(boid, boids, neighborRadius);
							forceSeparation += separation(boid, boids, avoidBoids);
							forceAlignment += alignment(boid, boids, neighborRadius);
							count++;
						}
					}
				}
			}
		}

		if (count > 0) {
			forceCohesion /= count;
			forceSeparation /= count;
			forceAlignment /= count;
		}

		glm::vec3 forceSeparationObstacles = separationObstacles(boid, obstacles, avoidObstacles);
		glm::vec3 forceSeparationBuildings = separationBuildings(boid, buildings) * 0.5f;
		glm::vec3 forceAttraction = attraction(boid);

		boid.acceleration = forceAlignment + forceCohesion + forceSeparation
			+ forceSeparationObstacles * 2 + forceSeparationBuildings * 5
			+ forceAttraction;
		boid.velocity += boid.acceleration * deltaTime * 0.8f;
		boid.velocity = limitSpeed(boid.velocity, maxSpeed);
		boid.position += boid.velocity * deltaTime;

		checkPosition(boid, minBoundary, maxBoundary);
		updateKDOPBoid(boid);
	}
}

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

void drawObjectColor(Core::RenderContext& context, glm::mat4 modelMatrix, glm::vec3 color) {

	glUseProgram(program);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(program, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(program, "color"), color.x, color.y, color.z);
	glUniform3f(glGetUniformLocation(program, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	Core::DrawContext(context);

}

void drawObjectTexture(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint texture) {

	glUseProgram(programTex);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(programTex, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(programTex, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(programTex, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	Core::SetActiveTexture(texture, "colorTexture", programTex, 0);
	Core::DrawContext(context);

}

void drawObjectNormal(Core::RenderContext& context, glm::mat4 modelMatrix, GLuint textureID, GLuint normalmapId, GLuint Program) {
	glUseProgram(Program);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(Program, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(Program, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(Program, "lightPos"), lightPos.x, lightPos.y, lightPos.z);

	Core::SetActiveTexture(textureID, "colorTexture", Program, 0);
	Core::SetActiveTexture(normalmapId, "normalSampler", Program, 1);
	Core::DrawContext(context);

}

void drawObjectBoid(Core::RenderContext& context, glm::mat4 modelMatrix, const Boid& boid) {
	glm::vec3 color = boid.color;
	glUseProgram(program);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(program, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(program, "color"), color.x, color.y, color.z);
	glUniform3f(glGetUniformLocation(program, "lightPos"), lightPos.x, lightPos.y, lightPos.z);
	Core::DrawContext(context);

}

void drawBoids() {
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

		drawObjectTexture(coneContext, modelMatrix, texture::dove);
	}
}

void drawObstacles() {
	for (const auto& obstacle : obstacles) {
		// Najpierw translacja
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), obstacle.position);

		// Następnie skalowanie
		modelMatrix = glm::scale(modelMatrix, glm::vec3(obstacle.size, obstacle.size, obstacle.size));

		// Rysowanie obiektu
		if (normalMappingEnabled) {
			drawObjectNormal(sphereContext, modelMatrix, texture::earth, texture::earthNormal, programNormal);
		}
		else {
			drawObjectTexture(sphereContext, modelMatrix, texture::earth);
		}
	}
}

void drawBuildings() {
	for (const auto& building : buildings) {
		// Najpierw translacja
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), building.position);

		// Następnie skalowanie
		modelMatrix = glm::scale(modelMatrix, glm::vec3(building.size.x, building.size.y, building.size.z));

		if (normalMappingEnabled) {
			drawObjectNormal(buildingContext, modelMatrix, texture::steel, texture::steelNormal, programNormal);
		}
		else {
			drawObjectTexture(buildingContext, modelMatrix, texture::steel);
		}
	}
}

void makescene() {

	drawBuildings();
	renderTerrain(vertices, indices);

}

void renderScene(GLFWwindow* window)
{
	glClearColor(0.0f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glm::mat4 view = createCameraMatrix();
	glm::mat4 projection = createPerspectiveMatrix();
	drawSkybox(view, projection);

	glUseProgram(program);

	double currentTime = glfwGetTime();
	double time = currentTime - lastTime;
	updateBoids(time, neighborRadius, avoidBoids);
	lastTime = currentTime;
	drawBoids();
	makescene();
	drawObstacles();

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
	glViewport(0, 0, width, height);
}

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
	
	glEnable(GL_DEPTH_TEST);
	program = shaderLoader.CreateProgram("shaders/shader.vert", "shaders/shader.frag");
	programTex = shaderLoader.CreateProgram("shaders/shader_tex.vert", "shaders/shader_tex.frag");
	shadowShaderProgram = shaderLoader.CreateProgram("shaders/shader_shadow.vert", "shaders/shader_shadow.frag");
	programNormal = shaderLoader.CreateProgram("shaders/shader_normal.vert", "shaders/shader_normal.frag");
	loadModelToContext("./models/dove.obj", coneContext);
	loadModelToContext("./models/sphere.obj", sphereContext);
	loadModelToContext("./models/cuboid.obj", buildingContext);
	initializeBoids(amountOfBoids, glm::vec3(0.0, 1.0, 0.3));
	initializeBoids(amountOfBoids, glm::vec3(0.0, 0.0, 1.0));

	texture::earth = Core::LoadTexture("./textures/earth.png");
	texture::earthNormal = Core::LoadTexture("./textures/earth_normalmap.png");

	texture::steel = Core::LoadTexture("./textures/Steel_S.jpg");
	texture::steelNormal = Core::LoadTexture("./textures/Steel_N.jpg");

	texture::dove = Core::LoadTexture("./textures/DOVE.jpg");

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
	if (glfwGetKey(window, GLFW_KEY_K) == GLFW_PRESS) {
		Obstacle obstacle;
		obstacle.position = cameraPos + cameraDir; // Pozycja przeszkody
		obstacle.size = 0.1f;          // Przykładowy rozmiar
		obstacle.vertices = sphereContext.getVertices();
		obstacle.indices = sphereContext.getIndices();
		obstacle.setKDOP();
		obstacles.push_back(obstacle);
	}
	if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) {
		attract = (attract + 1) % 2;
	}
	if (glfwGetKey(window, GLFW_KEY_Y) == GLFW_PRESS) {
		shadowMappingEnabled = !shadowMappingEnabled;  // Zmiana stanu shadow mappingu
		std::cout << "Shadow Mapping " << (shadowMappingEnabled ? "Enabled" : "Disabled") << std::endl;
	}
}

void renderLoop(GLFWwindow* window) {
	InitImGui(window);

	generateTerrain();
	setupBuildings();
	createTerrainMesh(vertices, indices);

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