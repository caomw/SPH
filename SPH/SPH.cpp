#include "SPH.h"
#include "GLWidget3D.h"
#include <GL/GLU.h>

#define SQR(x)	((x) * (x))

SPH::SPH(float radius_of_particle, float radius_of_smooth, float pressure_factor, float viscosity_factor, float density_0, float dumping_factor, float container_width, float container_depth, float container_height, float deltaT) {
	this->radius_of_particle = radius_of_particle;
	this->radius_of_smooth = radius_of_smooth;
	this->mass_of_particle = 10;//density_0 * container_width * container_depth * container_height / num_particles;
	this->pressure_factor = pressure_factor;
	this->viscosity_factor = viscosity_factor;
	this->density_0 = density_0;
	this->dumping_factor = dumping_factor;
	this->container_width = container_width;
	this->container_depth = container_depth;
	this->container_height = container_height;
	this->deltaT = deltaT;

	// randomly generate particles
	for (float x = container_width * 0.3; x <= container_width * 0.5 - radius_of_particle; x += radius_of_particle * 2) {
		for (float y = -container_depth * 0.5 + radius_of_particle; y <= container_depth * 0.5 - radius_of_particle; y += radius_of_particle * 2) {
			for (float z = radius_of_particle; z < container_height * 0.5; z += radius_of_particle * 2) {
				Particle particle(QVector3D(x, y, z));
				int r = random() * 255;
				int g = random() * 255;
				int b = random() * 255;
				particle.color = QColor(r, g, b);
				particles.push_back(particle);
			}
		}
	}
}

void SPH::update() {
	collisionDetection();
	updateDensity();
	updateForce();
	updateVelocityAndPosition();
}

void SPH::updateDensity() {
	for (int i = 0; i < particles.size(); ++i) {
		particles[i].density = 0.0f;

		for (int k = 0; k < particles[i].neighbors.size(); ++k) {
			int j = particles[i].neighbors[k];
			QVector3D r = particles[i].position - particles[j].position;
			particles[i].density += mass_of_particle * W_poly6(r.length(), radius_of_smooth);
		}
	}
}

void SPH::updateForce() {
	for (int i = 0; i < particles.size(); ++i) {
		QVector3D F_pressure;
		QVector3D F_viscosity;

		for (int k = 0; k < particles[i].neighbors.size(); ++k) {
			int j = particles[i].neighbors[k];
			if (j == i) continue;
			QVector3D r = particles[i].position - particles[j].position;
			F_pressure += -particles[i].density * mass_of_particle * (pressure(particles[i].density) / SQR(particles[i].density) + pressure(particles[j].density) / SQR(particles[j].density)) * dW_spiky(r, radius_of_smooth);
			F_viscosity += viscosity_factor * mass_of_particle * (particles[j].velocity - particles[i].velocity) / particles[j].density * ddW_viscosity(r.length(), radius_of_smooth);
		}

		//particles[i].force = F_pressure + F_viscosity + particles[i].density * QVector3D(0, 0, -98000000);
		particles[i].force = F_pressure + /*F_viscosity + */ QVector3D(0, 0, -9800);
	}
}

void SPH::updateVelocityAndPosition() {
	for (int i = 0; i < particles.size(); ++i) {
		QVector3D new_velocity = particles[i].velocity + particles[i].force / particles[i].density * deltaT;
		particles[i].position += (particles[i].velocity + new_velocity) * 0.5 * deltaT;
		particles[i].velocity = new_velocity;
	}
}

void SPH::collisionDetection() {
	for (int i = 0; i < particles.size(); ++i) {
		// check z coordinates
		if (particles[i].position.z() < radius_of_particle) {
			particles[i].position.setZ(radius_of_particle);
			particles[i].velocity.setZ(-particles[i].velocity.z() * dumping_factor);
		} else if (particles[i].position.z() > container_height - radius_of_particle) {
			particles[i].position.setZ(container_height - radius_of_particle);
			particles[i].velocity.setZ(-particles[i].velocity.z() * dumping_factor);
		}

		// check x coordinates
		if (particles[i].position.x() < -container_width * 0.5 + radius_of_particle) {
			particles[i].position.setX(-container_width * 0.5 + radius_of_particle);
			particles[i].velocity.setX(-particles[i].velocity.x() * dumping_factor);
		} else if (particles[i].position.x() > container_width * 0.5 - radius_of_particle) {
			particles[i].position.setX(container_width * 0.5 - radius_of_particle);
			particles[i].velocity.setX(-particles[i].velocity.x() * dumping_factor);
		}

		// check y coordinates
		if (particles[i].position.y() < -container_depth * 0.5 + radius_of_particle) {
			particles[i].position.setY(-container_depth * 0.5 + radius_of_particle);
			particles[i].velocity.setY(-particles[i].velocity.y() * dumping_factor);
		} else if (particles[i].position.y() > container_depth * 0.5 - radius_of_particle) {
			particles[i].position.setY(container_depth * 0.5 - radius_of_particle);
			particles[i].velocity.setY(-particles[i].velocity.y() * dumping_factor);
		}

		// check with others
		particles[i].neighbors.clear();
		for (int j = 0; j < particles.size(); ++j) {
			float r = (particles[i].position - particles[j].position).length();
			if (r < radius_of_smooth) {
				particles[i].neighbors.push_back(j);
		
				if (i != j && r < radius_of_particle * 2) {
					QVector3D n = particles[i].position - particles[j].position;
					n.normalize();
					particles[i].position = particles[j].position + n * radius_of_particle * 2;
					particles[i].velocity -= n * QVector3D::dotProduct(particles[i].velocity, n) * (1 + dumping_factor);
				}
			}
		}
	}
}


float SPH::W_poly6(float r, float h) {
	return 315.0f / 64 / M_PI / powf(h, 9) * powf(SQR(h) - SQR(r), 3);
}

QVector3D SPH::dW_spiky(const QVector3D& r, float h) {
	return -45 * SQR(h - r.length()) / M_PI / powf(h, 6) / r.length() * r;
}

float SPH::ddW_viscosity(float r, float h) {
	return 15.0f / 2 / M_PI / h / h / h * (-3 * r / 2 / h / h / h + 3 * h / 2 / r / r / r + 1 / r);
}

float SPH::pressure(float density) {
	return pressure_factor * pressure_factor * (density - density_0);
}

void SPH::draw() {
	for (int i = 0; i < particles.size(); ++i) {
		drawSphere(particles[i].position.x(), particles[i].position.y(), particles[i].position.z(), radius_of_particle, particles[i].color);
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