#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <algorithm>
#include <random>
#include <omp.h>
#include <chrono>
#include <iostream>

const unsigned WIN_W = 900;
const unsigned WIN_H = 700;

// fluid paramaters
const int   totalParticles = 5000; 
const float particleRadius = 5.5f;
const float particleMass   = 19.0f;
const float restingDensity = 1000.0f;
const float gasConstant    = 2000.0f;
const float viscosity      = 250.0f;
const float timestep       = 1.0f / 350.0f;
const float boundDamping   = -0.5f;
const float maxVelocity    = 600.0f;
const float PI_F = 3.14159265358979323846f;
float gravityForce = 800.0f;
float radSq, radFifth, radEighth;


struct Particle 
{
    float pressure;
    float density;
    sf::Vector2f currentForce, velocity, position;
};

std::vector<Particle> particles;

//use mouse to move fluid
struct MouseState 
{
    sf::Vector2f position, prevPosition;
    bool active = false;
    bool obstacleMode = false;
};
MouseState mouse;

//SPATIAL HASHING YEAH!!!
//particleCellIndices is where are particles actually live (their cell)
//particleIndices acts as a pointers to those cells (we sort this array to easily check neighbors)
//cellOffsets bookmarks our particleIndices so we know where each cell begins and ends
std::vector<uint32_t> particleCellIndices;
std::vector<uint32_t> particleIndices;
std::vector<uint32_t> cellOffsets;


//comparer for particle cell index from hash
struct ParticleIndexComparer 
{
    bool operator()(uint32_t x, uint32_t y) const 
    {
        return particleCellIndices[x] < particleCellIndices[y];
    }
};

ParticleIndexComparer comparer;

uint32_t HashCell(int cellX, int cellY) 
{
    const uint32_t primeX = 73856093;
    const uint32_t primeY = 19349663;
    
    uint32_t hash = ((uint32_t)cellX * primeX) ^ ((uint32_t)cellY * primeY);
    return hash % totalParticles;
}

void HashParticles() 
{
    for (int i = 0; i < totalParticles; i++) 
    {
        cellOffsets[i] = -1;

        int cellX = (int)(particles[i].position.x / particleRadius);
        int cellY = (int)(particles[i].position.y / particleRadius);
        
        particleCellIndices[i] = HashCell(cellX, cellY);
        //reset array
        particleIndices[i] = i; 
    }
}
//particles with the same cell are sorted and now contiguous in memory
void SortParticles() 
{
    std::sort(particleIndices.begin(), particleIndices.end(), comparer);
}

//lets us jump to cells needed to check neighbor particles
void CalculateCellOffsets() 
{
    for (int i = 0; i < totalParticles; i++) 
    {
        uint32_t particleIndex = particleIndices[i];
        uint32_t cellIndex = particleCellIndices[particleIndex];

        if (i == 0 || cellIndex != particleCellIndices[particleIndices[i - 1]]) 
        {
            cellOffsets[cellIndex] = i;
        }
    }
}

//kernels for our physics
//shoutout muller
//density
float StdKernel(float distSq) 
{
    if (distSq >= radSq) return 0.0f;
    float x = radSq - distSq;
    return (4.0f / (PI_F * radEighth)) * x * x * x;
}
//pressure
float SpikyKernelGradientMult(float distance) 
{
    float hr = particleRadius - distance;
    return (-10.0f / (PI_F * radFifth)) * hr * hr;
}
//viscocity
float ViscosityLaplacian(float distance) 
{
    float hr = particleRadius - distance;
    return (40.0f / (PI_F * radFifth)) * hr;
}
// missing surface tension but its ok

void ComputeDensityPressure() 
{
    //multithreaded
    #pragma omp parallel for
    for (int i = 0; i < totalParticles; i++) 
    {
        uint32_t particleIndex = particleIndices[i];
        
        int cellX = (int)(particles[particleIndex].position.x / particleRadius);
        int cellY = (int)(particles[particleIndex].position.y / particleRadius);
        float sum = 0.0f;

        //we do a 3x3 search area for neighbor calculations
        for (int x = -1; x <= 1; x++) 
        {
            for (int y = -1; y <= 1; y++) 
            {
                uint32_t targetCell = HashCell(cellX + x, cellY + y);
                uint32_t iterator = cellOffsets[targetCell];

                while (iterator != -1 && iterator < totalParticles) 
                {
                    uint32_t neighborIdx = particleIndices[iterator];
                    if (particleCellIndices[neighborIdx] != targetCell) break;

                    sf::Vector2f diff = particles[neighborIdx].position - particles[particleIndex].position;
                    float distSq = diff.x * diff.x + diff.y * diff.y;

                    if (distSq < radSq) 
                    {
                        sum += particleMass * StdKernel(distSq);
                    }
                    iterator++;
                }
            }
        }
        particles[particleIndex].density = sum + 0.000001f;
        particles[particleIndex].pressure = gasConstant * (particles[particleIndex].density - restingDensity);
    }
}

void ComputeForces() 
{
    //NOTICE: we recalculate the 9cell hash lookups
    //however, this is better than storing the entire list of neighbors in memory
    #pragma omp parallel for
    for (int i = 0; i < totalParticles; i++) 
    {
        uint32_t particleIndex = particleIndices[i];
        
        int cellX = (int)(particles[particleIndex].position.x / particleRadius);
        int cellY = (int)(particles[particleIndex].position.y / particleRadius);

        sf::Vector2f pressureForce(0.0f, 0.0f);
        sf::Vector2f viscosityForce(0.0f, 0.0f);

        //we do a 3x3 search area for neighbor calculations
        for (int x = -1; x <= 1; x++) 
        {
            for (int y = -1; y <= 1; y++) 
            {
                uint32_t targetCell = HashCell(cellX + x, cellY + y);
                uint32_t iterator = cellOffsets[targetCell];

                while (iterator != -1 && iterator < totalParticles) 
                {
                    uint32_t neighborIdx = particleIndices[iterator];
                    if (particleCellIndices[neighborIdx] != targetCell) break;

                    if (particleIndex != neighborIdx) 
                    {
                        sf::Vector2f diff = particles[neighborIdx].position - particles[particleIndex].position;
                        float distSq = diff.x * diff.x + diff.y * diff.y;

                        if (distSq < radSq && distSq > 0.0001f) 
                        {
                            float distance = std::sqrt(distSq);
                            sf::Vector2f direction = diff / distance;

                            float pForceMag = -particleMass * (particles[particleIndex].pressure + particles[neighborIdx].pressure)
                             / (2.0f * particles[neighborIdx].density) * SpikyKernelGradientMult(distance);
                            pressureForce += direction * pForceMag;

                            float vForceMag = viscosity * particleMass / particles[neighborIdx].density * ViscosityLaplacian(distance);
                            viscosityForce += (particles[neighborIdx].velocity - particles[particleIndex].velocity) * vForceMag;
                        }
                    }
                    iterator++;
                }
            }
        }

        sf::Vector2f gravity(0.0f, gravityForce * particles[particleIndex].density);
        particles[particleIndex].currentForce = pressureForce + viscosityForce + gravity;
        if (mouse.active && !mouse.obstacleMode) 
        {
            sf::Vector2f mouseVel = (mouse.position - mouse.prevPosition) * 30.0f;
            sf::Vector2f mouseDiff = particles[particleIndex].position - mouse.position;
            float mouseDistSq = mouseDiff.x * mouseDiff.x + mouseDiff.y * mouseDiff.y;
            
            if (mouseDistSq < 6400.0f) 
            {
                float f = 1.0f - std::sqrt(mouseDistSq) / 80.0f;
                particles[particleIndex].velocity += mouseVel * f * 0.1f;
            }
        }
    }
    if (mouse.active) mouse.prevPosition = mouse.position;
}

void Integrate() 
{
    #pragma omp parallel for
    for (int i = 0; i < totalParticles; i++) 
    {
        Particle& particle = particles[i];

        if (particle.density > 0.0001f) 
        {
            particle.velocity += (particle.currentForce / particle.density) * timestep;
        }

        //scale particle to maxVelocity vector if big
        float speedSq = particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y;
        if (speedSq > maxVelocity * maxVelocity) 
        {
            float speed = maxVelocity / std::sqrt(speedSq);
            particle.velocity *= speed;
        }

        particle.position += particle.velocity * timestep;
        // mouse force
        if (mouse.active && mouse.obstacleMode) 
        {
            const float mouseRadius = 60.0f;
            sf::Vector2f diff = particle.position - mouse.position;
            float distSq = diff.x * diff.x + diff.y * diff.y;
            
            //if particle is inside the mouse circle
            if (distSq < mouseRadius * mouseRadius && distSq > 0.0001f) 
            {
                float dist = std::sqrt(distSq);
                sf::Vector2f normal = diff / dist;
  
                particle.position = mouse.position + normal * mouseRadius;
                
                float velocityAlongNormal = particle.velocity.x * normal.x + particle.velocity.y * normal.y;
                if (velocityAlongNormal < 0.0f)
                {
                    particle.velocity -= normal * velocityAlongNormal * 1.5f; 
                }
            }
        }

        // box boundaries
        if (particle.position.x - particleRadius < 0.0f) 
        { 
            particle.velocity.x *= boundDamping; particle.position.x = particleRadius; 
        } else if (particle.position.x + particleRadius > WIN_W) 
        { 
            particle.velocity.x *= boundDamping; particle.position.x = WIN_W - particleRadius; 
        }

        if (particle.position.y - particleRadius < 0.0f) 
        { 
            particle.velocity.y *= boundDamping; particle.position.y = particleRadius; 
        } else if (particle.position.y + particleRadius > WIN_H) 
        { 
            particle.velocity.y *= boundDamping; particle.position.y = WIN_H - particleRadius; 
        }
    }
}

void FixedUpdate() 
{
    //setup spatial grid
    HashParticles();
    SortParticles();
    CalculateCellOffsets();

    ComputeDensityPressure();
    ComputeForces();
    Integrate();
}

void SpawnParticlesInBox() 
{
    particles.clear();
    particles.resize(totalParticles);
    gravityForce = 800.0f;
    //I dont want them uniform so spread randomly in spawn
    std::mt19937 rng(1337);
    std::uniform_real_distribution<float> jitter(-0.25f, 0.25f);

    float spacing = particleRadius * 0.90f;
    int cols = (int)std::ceil(std::sqrt((float)totalParticles * ((float)WIN_W / (float)WIN_H)));
    float startX = WIN_W * 0.25f;
    float startY = WIN_H * 0.10f;

    int i = 0;
    for (int row = 0; i < totalParticles; ++row) 
    {
        for (int col = 0; col < cols && i < totalParticles; ++col) 
        {
            particles[i].position.x = startX + col * spacing + jitter(rng);
            particles[i].position.y = startY + row * spacing;
            ++i;
        }
    }
}

int main() 
{
    radSq = particleRadius * particleRadius;
    radFifth = std::pow(particleRadius, 5);
    radEighth = std::pow(particleRadius, 8);

    particleIndices.resize(totalParticles);
    particleCellIndices.resize(totalParticles);
    cellOffsets.resize(totalParticles);

    sf::ContextSettings settings;
    settings.antialiasingLevel = 8;
    sf::RenderWindow window(sf::VideoMode(WIN_W, WIN_H), "awesome fluid sim :)", sf::Style::Default, settings);
    window.setFramerateLimit(60);
    SpawnParticlesInBox();

    //generate blurred particle texture, used for fluid shader
    //64x64 texture used at every particle position
    bool useShader = true;
    const int texSize = 64;
    sf::Image blurImg;
    blurImg.create(texSize, texSize, sf::Color::Transparent);
    for (int y = 0; y < texSize; ++y) 
    {
        for (int x = 0; x < texSize; ++x) 
        {
            float dx = x - (texSize / 2.0f);
            float dy = y - (texSize / 2.0f);
            float distSq = dx * dx + dy * dy;
            float radiusSq = (texSize / 2.0f) * (texSize / 2.0f);
            //only color pixels inside circle
            if (distSq < radiusSq) 
            {
                //gaussian blending
                //change param for softer or sharper circle
                float a = std::exp(-1.8f * distSq / radiusSq);
                blurImg.setPixel(x, y, sf::Color(255, 255, 255, (sf::Uint8)(a * 255)));
            }
        }
    }

    //uploading our texture to gpu
    sf::Texture particleTex;
    particleTex.loadFromImage(blurImg);
    
    sf::Sprite particleSprite(particleTex);
    particleSprite.setOrigin(texSize / 2.0f, texSize / 2.0f);
    //scaling texture down
    particleSprite.setScale(30.0f / texSize, 30.0f / texSize); 
    //render our particle blobs first onto this texture to feed into shader
    sf::RenderTexture fluidCanvas;
    fluidCanvas.create(WIN_W, WIN_H);

    //compile shaders
    sf::Shader shaderBasic;
    if (!shaderBasic.loadFromFile("../Resources/basic_shader.frag", sf::Shader::Fragment)) 
    {
        std::cerr << "Failed to load basic shader" << std::endl;
    }

    sf::Shader shaderAdvanced;
    if (!shaderAdvanced.loadFromFile("../Resources/advanced_shader.frag", sf::Shader::Fragment)) 
    {
        std::cerr << "Failed to load advanced shader" << std::endl;
    }
    int renderMode = 0;       // 0 = Dots, 1 = Basic Shader, 2 = Advanced Shader
    int currentBgIndex = 0;

    //setting up backgrounds
    std::vector<sf::Texture> bgTextures;

    sf::Image solidImg;
    solidImg.create(WIN_W, WIN_H, sf::Color(15, 25, 45));
    sf::Texture solidTex;
    solidTex.loadFromImage(solidImg);
    bgTextures.push_back(solidTex);

    sf::Texture bg1;
    if (bg1.loadFromFile("../Resources/backgrounds/brick.jpg")) bgTextures.push_back(bg1);

    sf::Texture bg4;
    if (bg4.loadFromFile("../Resources/backgrounds/coral.jpg")) bgTextures.push_back(bg4);

    sf::Texture bg3;
    if (bg3.loadFromFile("../Resources/backgrounds/check.png")) bgTextures.push_back(bg3);

    sf::Texture bg2;
    if (bg2.loadFromFile("../Resources/backgrounds/forest.jpg")) bgTextures.push_back(bg2);

    shaderBasic.setUniform("texture", sf::Shader::CurrentTexture);
    shaderBasic.setUniform("threshold", 0.05f);
    shaderBasic.setUniform("resolution", sf::Vector2f((float)WIN_W, (float)WIN_H));

    shaderAdvanced.setUniform("texture", sf::Shader::CurrentTexture);
    shaderAdvanced.setUniform("threshold", 0.05f);
    shaderAdvanced.setUniform("resolution", sf::Vector2f((float)WIN_W, (float)WIN_H));

    //additive blending
    sf::BlendMode fluidBlend(
        sf::BlendMode::SrcAlpha, sf::BlendMode::One, sf::BlendMode::Add, // rgb
        sf::BlendMode::One,      sf::BlendMode::One, sf::BlendMode::Add  // alpha, dense regions = +alpha
    );

    
    

    sf::CircleShape dot(3.0f);
    dot.setOrigin(3.0f, 3.0f);
    dot.setOutlineThickness(0.0f);

    sf::Clock clock;
    float acc = 0.0f;
    float totalPhysicsTime = 0.0f;
    int benchmarkFrameCount = 0;
    bool slowMode = false;

    while (window.isOpen()) 
    {
        sf::Event ev;
        while (window.pollEvent(ev)) 
        {
            if (ev.type == sf::Event::Closed) window.close();

            else if (ev.type == sf::Event::KeyPressed) 
            {
                if (ev.key.code == sf::Keyboard::Escape) window.close();
                else if (ev.key.code == sf::Keyboard::R) SpawnParticlesInBox();
                else if (ev.key.code == sf::Keyboard::S) slowMode = !slowMode;
                else if (ev.key.code == sf::Keyboard::Up) gravityForce = -800.0f;
                else if (ev.key.code == sf::Keyboard::Down) gravityForce = 800.0f;
                else if (ev.key.code == sf::Keyboard::M) mouse.obstacleMode = !mouse.obstacleMode;
                else if (ev.key.code == sf::Keyboard::Num1) renderMode = 0;
                else if (ev.key.code == sf::Keyboard::Num2) renderMode = 1;
                else if (ev.key.code == sf::Keyboard::Num3) renderMode = 2;
                else if (ev.key.code == sf::Keyboard::Left) 
                {
                    currentBgIndex--;
                    if (currentBgIndex < 0) currentBgIndex = bgTextures.size() - 1; 
                }
                else if (ev.key.code == sf::Keyboard::Right) 
                {
                    currentBgIndex++;
                    if (currentBgIndex >= bgTextures.size()) currentBgIndex = 0;
                }
            }
            
            else if (ev.type == sf::Event::MouseButtonPressed && ev.mouseButton.button == sf::Mouse::Left) 
            {
                mouse.position.x = mouse.prevPosition.x = (float)ev.mouseButton.x;
                mouse.position.y = mouse.prevPosition.y = (float)ev.mouseButton.y;
                mouse.active = true;
            }
            else if (ev.type == sf::Event::MouseButtonReleased && ev.mouseButton.button == sf::Mouse::Left) 
            {
                mouse.active = false;
            }
            else if (ev.type == sf::Event::MouseMoved && mouse.active) 
            {
                mouse.position.x = (float)ev.mouseMove.x;
                mouse.position.y = (float)ev.mouseMove.y;
            }
        }

        const float frameDt = std::min(1.0f / 30.0f, clock.restart().asSeconds());
        acc += frameDt;
        int steps = 0;
        int stepLimit;
        if(slowMode)
        {
            stepLimit = 1;
        }
        else
        {
            stepLimit = 8;
        }
        
        auto startTime = std::chrono::high_resolution_clock::now();
        while (acc >= timestep && steps < stepLimit) 
        {
            FixedUpdate();
            acc -= timestep;
            ++steps;
        }
        auto endTime = std::chrono::high_resolution_clock::now();

        if (steps > 0) 
        {
            std::chrono::duration<float, std::milli> elapsed = endTime - startTime;
            totalPhysicsTime += elapsed.count();
            benchmarkFrameCount += steps;
        }
        
        if (benchmarkFrameCount >= 500) 
        {
            float avgTime = totalPhysicsTime / (float)benchmarkFrameCount;
            std::cout << "Particles: " << totalParticles << " | Avg Time per Physics Step: " << avgTime << " ms" << std::endl;
            
            totalPhysicsTime = 0.0f;
            benchmarkFrameCount = 0;
        }

        //RENDERING
        window.clear();
        //first draw current background
        sf::Sprite bgSprite(bgTextures[currentBgIndex]);
        float scaleX = (float)WIN_W / bgTextures[currentBgIndex].getSize().x;
        float scaleY = (float)WIN_H / bgTextures[currentBgIndex].getSize().y;
        bgSprite.setScale(scaleX, scaleY);
        window.draw(bgSprite);

        
        
        if (renderMode == 0) //dots
        {
            for (const auto& particle : particles) 
            {
                float speed = std::sqrt(particle.velocity.x * particle.velocity.x + particle.velocity.y * particle.velocity.y);
                float t = std::min(1.0f, speed / 1000.0f);
                
                sf::Uint8 r = (sf::Uint8)(100 + t * 175);
                sf::Uint8 g = (sf::Uint8)(200 - t * 180);
                sf::Uint8 b = (sf::Uint8)(255 - t * 50);
                
                dot.setFillColor(sf::Color(r, g, b));
                dot.setPosition(particle.position.x, particle.position.y);
                window.draw(dot);
            }
        } 
        else //SHADERMODE 1 or 2
        {
            fluidCanvas.clear(sf::Color::Transparent);
            particleSprite.setColor(sf::Color(255, 255, 255, 15)); 
            
            for (const auto& particle : particles) 
            {
                particleSprite.setPosition(particle.position.x, particle.position.y);
                fluidCanvas.draw(particleSprite, fluidBlend);
            }
            fluidCanvas.display();
            sf::Sprite fullscreenQuad(fluidCanvas.getTexture());
            
            if (renderMode == 1) 
            {
                window.draw(fullscreenQuad, &shaderBasic);
            } 
            else if (renderMode == 2) 
            {
                shaderAdvanced.setUniform("background", bgTextures[currentBgIndex]);
                window.draw(fullscreenQuad, &shaderAdvanced);
            }
        }

        //drawing mouse visual
        if (mouse.active) 
        {
            sf::CircleShape mouseCircle;
            if (mouse.obstacleMode) 
            {
                //solid mode
                mouseCircle.setRadius(60.0f);
                mouseCircle.setOrigin(60.0f, 60.0f);
                mouseCircle.setFillColor(sf::Color(235, 235, 235));
                mouseCircle.setOutlineThickness(2.0f);
                mouseCircle.setOutlineColor(sf::Color(200, 200, 200));
            } 
            else 
            {
                //pulling mode
                mouseCircle.setRadius(80.0f);
                mouseCircle.setOrigin(80.0f, 80.0f);
                mouseCircle.setFillColor(sf::Color(0, 50, 50, 100)); 
                mouseCircle.setOutlineThickness(1.0f);
                mouseCircle.setOutlineColor(sf::Color(100, 150, 255, 150));
            }
            
            mouseCircle.setPosition(mouse.position.x, mouse.position.y);
            window.draw(mouseCircle);
        }
        window.display();
    }
    return 0;
}
