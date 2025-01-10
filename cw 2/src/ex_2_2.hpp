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
glm::vec3 cameraDir =  glm::vec3(1.f, 0.f, 0.f);;
glm::vec3 spaceshipPos = glm::vec3(-4.f, 0, 0);
glm::vec3 spaceshipDir = glm::vec3(1.f, 0.f, 0.f);

float aspectRatio = 1.f;


float neighborRadius = 1.0f;
float sightAngle = 120.0f;
float avoidRadius = 0.4f;
float avoidRadius2 = 1.0f;
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

void initializeBoids(int numBoids, glm::vec3 color) {
	for (int i = 0; i < numBoids; i++) {
		Boid boid;
		boid.position = glm::vec3(rand() % 5/10.0f - 2.5f , rand() % 10 /10.0f, rand() % 10 / 10.0f);
		boid.velocity = glm::vec3((rand() % 20 - 10) / 10.0f, (rand() % 20 - 10) / 10.0f, (rand() % 20 - 10) / 10.0f);
		boid.acceleration = glm::vec3(0.0f);
		boid.color = color;
		boids.push_back(boid);
	}
}

bool checkCollision(const glm::vec3& min1, const glm::vec3& max1, const glm::vec3& min2, const glm::vec3& max2) {
	// Sprawdzamy, czy zakresy na każdej z osi się nie nakładają
	if (max1.x < min2.x || min1.x > max2.x) return false; // Oś X
	if (max1.y < min2.y || min1.y > max2.y) return false; // Oś Y
	if (max1.z < min2.z || min1.z > max2.z) return false; // Oś Z

	// Jeśli żadna z osi nie jest rozłączna, to kolizja występuje
	return true;
}

void checkCollisionBoid(Boid &boid) {
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

glm::vec3 separation(Boid& boid, const std::vector<Boid>& boids, float avoidRadius) {
	glm::vec3 avoid(0.0f);  // Siła unikania
	int count = 0;

	// Przeszukaj wszystkie boidy
	for (const auto& other : boids) {
		// Upewnijmy się, że boidy są w zasięgu widzenia i w obrębie unikania
		if (insight(boid, other.position, avoidRadius) && boid.color == other.color) {
			float distance = glm::distance(boid.position, other.position);
			if (distance < avoidRadius) {
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

glm::vec3 separationObstacles(Boid& boid, const std::vector<Obstacle>& obstacles, float avoidRadius2) {
	glm::vec3 avoid(0.0f);  // Siła unikania
	int count = 0;


	for (const auto& other : obstacles) {
		// Upewnijmy się, że boidy są w zasięgu widzenia i w obrębie unikania
		if (insight(boid, other.position, avoidRadius2)) {
			float distance = glm::distance(boid.position, other.position);
			if (distance < avoidRadius) {
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

void updateboxBoid(Boid &boid, glm::vec3 boxsize) {
	boid.minbox = boid.position - boxsize;
	boid.maxbox = boid.position + boxsize;
}


void updateBoids(float deltaTime, float neighborRadius, float avoidRadius) {
	
	for (auto& boid : boids) {
		glm::vec3 forceCohesion = cohesion(boid, boids, neighborRadius);
		glm::vec3 forceSeparation = separation(boid, boids,avoidRadius);
		glm::vec3 forceSeparationObstacles = separationObstacles(boid, obstacles, avoidRadius2);
		glm::vec3 forceSeparationBuildings = separationBuildings(boid, buildings);
		glm::vec3 forceAlignment = alignment(boid, boids, neighborRadius);
		glm::vec3 forceattraction = attraction(boid);

		boid.acceleration =  forceAlignment + forceCohesion  + forceSeparation + forceSeparationObstacles * 5 + forceSeparationBuildings * 10 + forceattraction ;
		boid.velocity += boid.acceleration * deltaTime;
		boid.velocity = limitSpeed(boid.velocity, maxSpeed);
		boid.position += boid.velocity * deltaTime;

		glm::vec3 minBoundary(-2.f, 0.0f, -2.0f);
		glm::vec3 maxBoundary(2.f, 2.0f, 2.0f);
		checkPosition(boid, minBoundary, maxBoundary);

		updateboxBoid(boid, boxsize);
		checkCollisionBoid(boid);
	}
}

glm::mat4 createPerspectiveMatrix()
{

	glm::mat4 perspectiveMatrix;
	float n = 0.05;
	float f = 20.;
	float a1 = glm::min(aspectRatio, 1.f);
	float a2 = glm::min(1 / aspectRatio, 1.f);
	perspectiveMatrix = glm::mat4({
		1,0.,0.,0.,
		0.,aspectRatio,0.,0.,
		0.,0.,(f + n) / (n - f),2 * f * n / (n - f),
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
	glUniform3f(glGetUniformLocation(program, "lightPos"), 0, 1, 0);
	Core::DrawContext(context);

}

void drawObjectBoid(Core::RenderContext& context, glm::mat4 modelMatrix,const Boid& boid) {
	glm::vec3 color = boid.color;
	glUseProgram(program);
	glm::mat4 viewProjectionMatrix = createPerspectiveMatrix() * createCameraMatrix();
	glm::mat4 transformation = viewProjectionMatrix * modelMatrix;
	glUniformMatrix4fv(glGetUniformLocation(program, "transformation"), 1, GL_FALSE, (float*)&transformation);
	glUniformMatrix4fv(glGetUniformLocation(program, "modelMatrix"), 1, GL_FALSE, (float*)&modelMatrix);
	glUniform3f(glGetUniformLocation(program, "color"), color.x, color.y, color.z);
	glUniform3f(glGetUniformLocation(program, "lightPos"), 0, 1, 0);
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
		glm::mat4 translationMatrix= glm::translate(glm::mat4(1.0f), boid.position);


		glm::mat4 modelMatrix = translationMatrix * rotationMatrix;

		drawObjectBoid(coneContext, modelMatrix, boid);
	}
}

void drawObstacles() {
	for (const auto& obstacle : obstacles) {
		// Najpierw translacja
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), obstacle.position);

		// Następnie skalowanie
		modelMatrix = glm::scale(modelMatrix, glm::vec3(obstacle.size, obstacle.size, obstacle.size));

		// Rysowanie obiektu
		drawObjectColor(sphereContext, modelMatrix, glm::vec3(1.0, 0.0, 0.0));
	}
}

void drawBuildings() {
	
	
	for (const auto& building : buildings) {
		// Najpierw translacja
		glm::mat4 modelMatrix = glm::translate(glm::mat4(1.0f), building.position);

		// Następnie skalowanie
		modelMatrix = glm::scale(modelMatrix, glm::vec3(building.size.x, building.size.y, building.size.z));

		// Rysowanie obiektu
		drawObjectColor(buildingContext, modelMatrix, glm::vec3(1.0, 1.0, 0.0));
	}
}

void drawSpaceship(const glm::mat4& cameraMatrix, glm::vec3 cameraDir, glm::vec3 cameraPos) {
	// Wektory boczny i górny statku, bazujące na kamerze
	glm::vec3 spaceshipSide = glm::normalize(glm::cross(cameraDir, glm::vec3(0.f, 1.f, 0.f)));
	glm::vec3 spaceshipUp = glm::normalize(glm::cross(spaceshipSide, cameraDir));

	// Macierz rotacji statku (zgodna z kamerą)
	glm::mat4 spaceshipRotationMatrix = glm::mat4({
		spaceshipSide.x, spaceshipSide.y, spaceshipSide.z, 0,
		spaceshipUp.x, spaceshipUp.y, spaceshipUp.z, 0,
		-cameraDir.x, -cameraDir.y, -cameraDir.z, 0,
		0.f, 0.f, 0.f, 1.f,
		});

	// Dodaj obrót o 180 stopni wokół osi Y, aby obrócić statek
	glm::mat4 rotation180 = glm::rotate(glm::mat4(1.0f), glm::radians(180.0f), glm::vec3(0.0f, 1.0f, 0.0f));

	// Ustaw pozycję statku (przesunięcie od pozycji kamery, np. 1.5 jednostki w przód)
	glm::vec3 spaceshipPosition = cameraPos + 1.5f * cameraDir;

	// Macierz modelu statku
	glm::mat4 spaceshipModelMatrix = glm::translate(spaceshipPosition) * spaceshipRotationMatrix * rotation180;

	// Narysuj statek
	drawObjectColor(shipContext, spaceshipModelMatrix, glm::vec3(0.3, 0.3, 0.5));
}


void makescene() {

	drawBuildings();
}

void renderScene(GLFWwindow* window)
{
	glClearColor(0.0f, 0.3f, 0.3f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(program);

	double currentTime = glfwGetTime();
	double time = currentTime - lastTime;
	updateBoids(time, neighborRadius, avoidRadius);
	lastTime = currentTime;
	drawBoids();
	makescene();
	drawObstacles();

	glm::mat4 cameraMatrix = createCameraMatrix();
	drawSpaceship(cameraMatrix, cameraDir, cameraPos);
	

	glUseProgram(0);

	glfwSwapBuffers(window);
}
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
	glViewport(0, 0, width, height);
}

struct Vertex {
	glm::vec3 position;
	glm::vec3 normal;
	glm::vec2 texCoords;
};

void loadModelToContext(std::string path, Core::RenderContext& context)
{
	Assimp::Importer import;
	const aiScene * scene = import.ReadFile(path, aiProcess_Triangulate | aiProcess_CalcTangentSpace);

	if (!scene || scene->mFlags & AI_SCENE_FLAGS_INCOMPLETE || !scene->mRootNode)
	{
		std::cout << "ERROR::ASSIMP::" << import.GetErrorString() << std::endl;
		return;
	}
	context.initFromAssimpMesh(scene->mMeshes[0]);
}

void addBuilding(glm::vec3 buildPos) {
	Building building;
	building.position = buildPos; // Pozycja przeszkody
	building.size = glm::vec3(0.5f,1.0f,1.0f);          // Przykładowy rozmiar
	building.minbox = building.position - buildingboxsize;
	building.maxbox = building.position + buildingboxsize;
	buildings.push_back(building);
}
void init(GLFWwindow* window)
{
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);

	glEnable(GL_DEPTH_TEST);
	program = shaderLoader.CreateProgram("shaders/shader_2_2.vert", "shaders/shader_2_2.frag");
	loadModelToContext("./models/cone.obj", coneContext);
	loadModelToContext("./models/sphere.obj", sphereContext);
	loadModelToContext("./models/cuboid.obj", buildingContext);
	loadModelToContext("./models/spaceship.obj", shipContext);
	initializeBoids(50, glm::vec3(0.0, 1.0, 0.3));
	initializeBoids(50, glm::vec3(0.0, 0.0, 1.0));
	addBuilding(glm::vec3(-1.0f, 0.5f, 0.0f));
	addBuilding(glm::vec3(-1.0f, 0.5f, 2.0f));
	addBuilding(glm::vec3(1.0f, 0.5f, 1.0f));
	addBuilding(glm::vec3(0.0f, 1.5f, 1.0f));

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
		cameraPos  += glm::vec3(0, 1, 0) * moveSpeed;
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
		obstacle.minbox = obstacle.position - obstacleboxsize;
		obstacle.maxbox = obstacle.position + obstacleboxsize;
		obstacles.push_back(obstacle);
	}
	if (glfwGetKey(window, GLFW_KEY_M) == GLFW_PRESS) {
		attract = (attract + 1) % 2;
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