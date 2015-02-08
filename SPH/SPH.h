#pragma once

#include <QVector3D>
#include <vector>
#include <QColor>

class Particle {
public:
	QVector3D position;
	QVector3D velocity;
	float density;
	QVector3D force;
	QColor color;
	std::vector<int> neighbors;
	std::vector<int> containers;

	Particle() {};
	Particle(const QVector3D& position) { this->position = position; }
};

class SPH {
private:
	float radius_of_particle;
	float radius_of_smooth;
	float mass_of_particle;
	float pressure_factor;
	float viscosity_factor;
	float density_0;
	std::vector<Particle> particles;
	float radius_of_container;
	float height_of_container;
	std::vector<Particle> containers;
	float deltaT;

public:
	SPH(int num_particles, float radius_of_particle, float radius_of_smooth, float pressure_factor, float viscosity_factor, float density_0, float radius_of_container, float height_of_container, float deltaT);

	void update();
	void updateDensity();
	void updateForce();
	void updateVelocityAndPosition();
	float W_poly6(float r, float h);
	QVector3D dW_spiky(const QVector3D& r, float h);
	float ddW_viscosity(float r, float h);
	float pressure(float density);
	void draw();
	void drawSphere(float x, float y, float z, float r, const QColor& color);
	float random();
};

