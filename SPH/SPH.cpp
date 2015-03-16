#include "SPH.h"
#include "GLWidget3D.h"
#include <GL/GLU.h>

#define SQR(x)	((x) * (x))

SPH::SPH(double radius_of_smooth, double mass_of_particle, double pressure_factor, double viscosity_factor, double density_0, double surface_coefficient, double surface_threshold, double dumping_factor, double container_width, double container_depth, double container_height, double deltaT) {
	this->radius_of_smooth = radius_of_smooth;
	//this->radius_of_smooth = pow(3.0 * 20.0 * mass_of_particle / 4.0 / M_PI / density_0, 1.0 / 3.0);
	this->mass_of_particle = mass_of_particle;
	this->pressure_factor = pressure_factor;
	this->viscosity_factor = viscosity_factor;
	this->density_0 = density_0;
	this->surface_coefficient = surface_coefficient;
	this->surface_threshold = surface_threshold;
	this->dumping_factor = dumping_factor;
	this->container_width = container_width;
	this->container_depth = container_depth;
	this->container_height = container_height;
	this->deltaT = deltaT;

	//float rho = 1000.0f;
	//float radius_of_particle = 0.5f / powf(rho / mass_of_particle, 1.0/3.0);
	float radius_of_particle = pow(0.02 / density_0, 1.0/3.0) * 0.5;


	// randomly generate particles
	for (float x = container_width * -0.5 + radius_of_particle; x <= container_width * -0.25; x += radius_of_particle * 2) {
		for (float y = -container_depth * 0.5 + radius_of_particle; y <= container_depth * 0.5 - radius_of_particle; y += radius_of_particle * 2) {
			for (float z = radius_of_particle; z < container_height - radius_of_particle; z += radius_of_particle * 2) {
				Particle particle(QVector3D(x, y, z));
				int r = 0;//random() * 255;
				int g = 0;//random() * 255;
				int b = 255;//random() * 255;
				particle.color = QColor(r, g, b);
				particles.push_back(particle);
			}
		}
	}
}

void SPH::update() {
	updateNeighbors();
	updateDensity();
	updateInternalForce();
	updateExternalForce();
	updateVelocityAndPosition();
	collisionDetection();
}

void SPH::updateDensity() {
	for (int i = 0; i < particles.size(); ++i) {
		particles[i].density = 0.0f;

		for (int k = 0; k < particles[i].neighbors.size(); ++k) {
			int j = particles[i].neighbors[k];
			QVector3D r = particles[i].position - particles[j].position;
			particles[i].density += mass_of_particle * W_poly6(r.length(), radius_of_smooth);
		}

		particles[i].pressure = pressure_factor * (particles[i].density - density_0);
	}
}

void SPH::updateInternalForce() {
	for (int i = 0; i < particles.size(); ++i) {
		QVector3D F_pressure;
		QVector3D F_viscosity;

		for (int k = 0; k < particles[i].neighbors.size(); ++k) {
			int j = particles[i].neighbors[k];
			if (j == i) continue;
			QVector3D r = particles[i].position - particles[j].position;
			F_pressure += -particles[i].density * mass_of_particle * (particles[i].pressure / SQR(particles[i].density) + particles[j].pressure / SQR(particles[j].density)) * dW_spiky(r, radius_of_smooth);
			F_viscosity += viscosity_factor * mass_of_particle * (particles[j].velocity - particles[i].velocity) / particles[j].density * ddW_viscosity(r.length(), radius_of_smooth);
		}

		particles[i].internalForce = F_pressure + F_viscosity;
	}
}

void SPH::updateExternalForce() {
	for (int i = 0; i < particles.size(); ++i) {
		QVector3D F_surface;
		
		QVector3D normal;
		for (int k = 0; k < particles[i].neighbors.size(); ++k) {
			int j = particles[i].neighbors[k];
			QVector3D r = particles[i].position - particles[j].position;

			normal += mass_of_particle / particles[j].density * dW_poly6(r, radius_of_smooth);
		}

		if (normal.length() > surface_threshold) {
			double c = 0.0;

			for (int k = 0; k < particles[i].neighbors.size(); ++k) {
				int j = particles[i].neighbors[k];
				QVector3D r = particles[i].position - particles[j].position;

				c += mass_of_particle / particles[j].density * ddW_poly6(r.length(), radius_of_smooth);
			}

			F_surface = -surface_coefficient * c * normal.normalized();
		}

		particles[i].externalForce = F_surface + particles[i].density * QVector3D(0, 0, -9.82);
	}
}

void SPH::updateVelocityAndPosition() {
	for (int i = 0; i < particles.size(); ++i) {
		QVector3D new_velocity = particles[i].velocity + (particles[i].internalForce + particles[i].externalForce) / particles[i].density * deltaT;
		particles[i].position += (particles[i].velocity + new_velocity) * 0.5 * deltaT;
		particles[i].velocity = new_velocity;
	}
}

void SPH::updateNeighbors() {
	for (int i = 0; i < particles.size(); ++i) {
		particles[i].neighbors.clear();

		for (int j = 0; j < particles.size(); ++j) {
			float r = (particles[i].position - particles[j].position).length();
			if (r < radius_of_smooth) {
				particles[i].neighbors.push_back(j);		
			}
		}
	}
}

void SPH::collisionDetection() {
	for (int i = 0; i < particles.size(); ++i) {
		// check z coordinates
		if (particles[i].position.z() < radius_of_particle) {
			respondCollision(particles[i].velocity, QVector3D(0, 0, 1), abs(particles[i].position.z() - radius_of_particle));
			particles[i].position.setZ(radius_of_particle);
		} else if (particles[i].position.z() > container_height - radius_of_particle) {
			respondCollision(particles[i].velocity, QVector3D(0, 0, -1), abs(particles[i].position.z() - container_height + radius_of_particle));
			particles[i].position.setZ(container_height - radius_of_particle);
		}

		// check x coordinates
		if (particles[i].position.x() < -container_width * 0.5 + radius_of_particle) {
			respondCollision(particles[i].velocity, QVector3D(1, 0, 0), abs(particles[i].position.x() + container_width * 0.5 - radius_of_particle));
			particles[i].position.setX(-container_width * 0.5 + radius_of_particle);
		} else if (particles[i].position.x() > container_width * 0.5 - radius_of_particle) {
			respondCollision(particles[i].velocity, QVector3D(-1, 0, 0), abs(particles[i].position.x() - container_width * 0.5 + radius_of_particle));
			particles[i].position.setX(container_width * 0.5 - radius_of_particle);
		}

		// check y coordinates
		if (particles[i].position.y() < -container_depth * 0.5 + radius_of_particle) {
			respondCollision(particles[i].velocity, QVector3D(0, 1, 0), abs(particles[i].position.y() + container_depth * 0.5 - radius_of_particle));
			particles[i].position.setY(-container_depth * 0.5 + radius_of_particle);
		} else if (particles[i].position.y() > container_depth * 0.5 - radius_of_particle) {
			respondCollision(particles[i].velocity, QVector3D(0, -1, 0), abs(particles[i].position.y() - container_depth * 0.5 + radius_of_particle));
			particles[i].position.setY(container_depth * 0.5 - radius_of_particle);
		}

		// check with others
		/*
		for (int k = 0; k < particles[i].neighbors.size(); ++k) {
			int j = particles[i].neighbors[k];
			if (i == j) continue;

			float r = (particles[i].position - particles[j].position).length();
			if (r < radius_of_particle) {
				QVector3D n = particles[i].position - particles[j].position;
				n.normalize();
				double d = abs(r - radius_of_particle);
				particles[i].position = particles[j].position + n * radius_of_particle * 2;
				particles[i].velocity -= n * QVector3D::dotProduct(particles[i].velocity, n) * (1 + dumping_factor * d / deltaT / particles[i].velocity.length());
			}
		}
		*/
	}
}

/**
 * 衝突に基づき、速度ベクトルを更新する。
 *
 * @param v [OUT]	現在の速度ベクトル ⇒ 更新後の速度ベクトルに更新される
 * @param n			衝突における法線ベクトル
 * @param d			めり込んだ距離
 */
void SPH::respondCollision(QVector3D& v, const QVector3D& n, double d) {
	v = v - (1.0 + dumping_factor * d / deltaT / v.length()) * QVector3D::dotProduct(v, n) * n;
}


double SPH::W_poly6(double r, double h) {
	return 315.0f / 64.0 / M_PI / pow(h, 9.0) * pow(SQR(h) - SQR(r), 3.0);
}

QVector3D SPH::dW_poly6(const QVector3D& r, double h) {
	return -945.0 / 32.0 / M_PI / pow(h, 9.0) * SQR(SQR(h) - r.lengthSquared()) * r;
}

double SPH::ddW_poly6(double r, double h) {
	return -945.0 / 32.0 / M_PI / pow(h, 9.0) * (SQR(h) - SQR(r)) * (3.0 * SQR(h) - 7.0 * SQR(r));
}

QVector3D SPH::dW_spiky(const QVector3D& r, double h) {
	double norm = r.length();
	if (fabs(norm) < 0.000001) {
		return -45.0 / M_PI / pow(h, 4.0) * r.normalized();
	} else {
		return -45.0 * SQR(h - norm) / M_PI / powf(h, 6) / norm * r;
	}
}

double SPH::ddW_viscosity(double r, double h) {
	return 45.0 / M_PI / pow(h, 6.0) * (h - r);
}

void SPH::draw() {
	for (int i = 0; i < particles.size(); ++i) {
		drawSphere(particles[i].position.x(), particles[i].position.y(), particles[i].position.z(), 0.02, particles[i].color);
	}
	for (int i = 0; i < containers.size(); ++i) {
		drawSphere(containers[i].position.x(), containers[i].position.y(), containers[i].position.z(), radius_of_particle, containers[i].color);
	}
}

void SPH::drawSphere(float x, float y, float z, float r, const QColor& color) {
	int slices = 16;
	int stacks = 8;

	glBegin(GL_QUADS);
	glColor3f(color.redF(), color.greenF(), color.blueF());
	for (int i = 0; i < slices; ++i) {
		float theta1 = M_PI * 2.0f / slices * i;
		float theta2 = M_PI * 2.0f / slices * (i + 1);

		for (int j = 0; j < stacks; ++j) {
			float phi1 = M_PI / stacks * j - M_PI * 0.5;
			float phi2 = M_PI / stacks * (j + 1) - M_PI * 0.5;

			QVector3D pt1 = QVector3D(cosf(theta1) * cosf(phi1), sinf(theta1) * cosf(phi1), sinf(phi1));
			QVector3D pt2 = QVector3D(cosf(theta2) * cosf(phi1), sinf(theta2) * cosf(phi1), sinf(phi1));
			QVector3D pt3 = QVector3D(cosf(theta2) * cosf(phi2), sinf(theta2) * cosf(phi2), sinf(phi2));
			QVector3D pt4 = QVector3D(cosf(theta1) * cosf(phi2), sinf(theta1) * cosf(phi2), sinf(phi2));

			glNormal3f(pt1.x(), pt1.y(), pt1.z());
			glVertex3f(x + pt1.x() * r, y + pt1.y() * r, z + pt1.z() * r);
			glNormal3f(pt2.x(), pt2.y(), pt2.z());
			glVertex3f(x + pt2.x() * r, y + pt2.y() * r, z + pt2.z() * r);
			glNormal3f(pt3.x(), pt3.y(), pt3.z());
			glVertex3f(x + pt3.x() * r, y + pt3.y() * r, z + pt3.z() * r);
			glNormal3f(pt4.x(), pt4.y(), pt4.z());
			glVertex3f(x + pt4.x() * r, y + pt4.y() * r, z + pt4.z() * r);
		}
	}
	glEnd();
}

float SPH::random() {
	return (float)(rand() % RAND_MAX) / RAND_MAX;
} 