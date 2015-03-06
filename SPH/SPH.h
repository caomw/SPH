#pragma once

#include <QVector3D>
#include <vector>
#include <QColor>

class Particle {
public:
	QVector3D position;
	QVector3D velocity;
	float density;
	float pressure;
	QVector3D internalForce;
	QVector3D externalForce;
	QColor color;
	std::vector<int> neighbors;
	std::vector<int> containers;

	Particle() {};
	Particle(const QVector3D& position) { this->position = position; }
};

class SPH {
private:
	double radius_of_particle;
	double radius_of_smooth;
	double mass_of_particle;
	double pressure_factor;
	double viscosity_factor;
	double density_0;
	double surface_coefficient;
	double surface_threshold;
	double dumping_factor;
	std::vector<Particle> particles;

	double container_width;
	double container_depth;
	double container_height;
	std::vector<Particle> containers;
	double deltaT;

public:
	SPH(double radius_of_smooth, double mass_of_particle, double pressure_factor, double viscosity_factor, double density_0, double surface_coefficient, double surface_threshold, double dumping_factor, double container_width, double container_depth, double container_height, double deltaT);

	void update();
	void updateDensity();
	void updateInternalForce();
	void updateExternalForce();
	void updateVelocityAndPosition();
	void updateNeighbors();
	void collisionDetection();
	double W_poly6(double r, double h);
	QVector3D dW_poly6(const QVector3D& r, double h);
	double ddW_poly6(double r, double h);
	QVector3D dW_spiky(const QVector3D& r, double h);
	double ddW_viscosity(double r, double h);
	void draw();
	void drawSphere(float x, float y, float z, float r, const QColor& color);
	float random();
};

