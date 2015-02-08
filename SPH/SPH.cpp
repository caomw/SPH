#include "SPH.h"
#include "GLWidget3D.h"
#include <GL/GLU.h>

SPH::SPH(int num_particles, float radius_of_particle, float radius_of_smooth, float pressure_factor, float viscosity_factor, float density_0, float radius_of_container, float height_of_container, float deltaT) {
	particles.resize(num_particles);
	this->radius_of_particle = radius_of_particle;
	this->radius_of_smooth = radius_of_smooth;
	this->mass_of_particle = radius_of_container * radius_of_container * M_PI * height_of_container / num_particles;
	this->pressure_factor = pressure_factor;
	this->viscosity_factor = viscosity_factor;
	this->density_0 = density_0;
	this->radius_of_container = radius_of_container;
	this->height_of_container = height_of_container;
	this->deltaT = deltaT;

	// randomly generate particles
	for (int i = 0; i < num_particles; ++i) {
		float x, y;
		while (true) {
			//x = radius_of_container * 2 * random() - radius_of_container;
			//y = radius_of_container * 2 * random() - radius_of_container;
			x = radius_of_container *  random() - radius_of_container * 0.5;
			y = radius_of_container *  random() - radius_of_container * 0.5;
			if (x * x + y * y < radius_of_container * radius_of_container) break;
		}
		float z = height_of_container * random();

		particles[i].position = QVector3D(x, y, z);

		int r = random() * 255;
		int g = random() * 255;
		int b = random() * 255;
		particles[i].color = QColor(r, g, b);
	}

	// build a container
	/*
	for (float x = -radius_of_container; x <= radius_of_container; x += radius_of_particle) {
		for (float y = -radius_of_container; y <= radius_of_container; y += radius_of_particle) {
			if (x * x + y * y > radius_of_container * radius_of_container) continue;

			Particle particle(QVector3D(x, y, 0));
			particle.color = QColor(128, 128, 128);
			containers.push_back(particle);
			//containers.push_back(Particle(QVector3D(x, y, height_of_container)));
		}
	}
	*/
	/*
	for (float z = radius_of_particle * 0.5; z < height_of_container; z += radius_of_particle * 0.5) {
		float d_theta = radius_of_container * 2 * M_PI / radius_of_particle * 2;
		for (float theta = 0.0f; theta < M_PI * 2; theta += d_theta) {
			float x = cosf(theta) * radius_of_container;
			float y = sinf(theta) * radius_of_container;
			containers.push_back(Particle(QVector3D(x, y, z)));
		}
	}
	*/
}

void SPH::update() {
	updateDensity();
	updateForce();
	updateVelocityAndPosition();
}

void SPH::updateDensity() {
	for (int i = 0; i < particles.size(); ++i) {
		particles[i].density = 0.0f;
		particles[i].neighbors.clear();

		for (int j = 0; j < particles.size(); ++j) {
			float r = (particles[i].position - particles[j].position).length();
			if (r < radius_of_smooth) {
				particles[i].neighbors.push_back(j);
				particles[i].density += mass_of_particle * W_poly6(r, radius_of_smooth);
			}
		}

		for (int j = 0; j < containers.size(); ++j) {
			float r = (particles[i].position - containers[j].position).length();
			if (r < radius_of_smooth) {
				particles[i].containers.push_back(j);
				float hoge = W_poly6(r, radius_of_smooth);
				float hoge2 = mass_of_particle * W_poly6(r, radius_of_smooth);
				particles[i].density += mass_of_particle * W_poly6(r, radius_of_smooth);
			}
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
			F_pressure += -mass_of_particle * (pressure(particles[i].density) + pressure(particles[j].density)) / 2 / particles[j].density * dW_spiky(r, radius_of_smooth);
			F_viscosity += viscosity_factor * mass_of_particle * (particles[j].velocity - particles[i].velocity) / particles[j].density * ddW_viscosity(r.length(), radius_of_smooth);
		}

		for (int k = 0; k < particles[i].containers.size(); ++k) {
			int j = particles[i].containers[k];
			QVector3D r = particles[i].position - containers[j].position;
			F_pressure += -mass_of_particle * (pressure(particles[i].density) + pressure(100)) / 2 / 100 * dW_spiky(r, radius_of_smooth);
			F_viscosity += viscosity_factor * mass_of_particle * - particles[i].velocity / 100 * ddW_viscosity(r.length(), radius_of_smooth);
		}

		//particles[i].force = F_pressure + F_viscosity + particles[i].density * QVector3D(0, 0, -98000000);
		particles[i].force = F_pressure * 0.001 + F_viscosity + QVector3D(0, 0, -98000000);
	}
}

void SPH::updateVelocityAndPosition() {
	for (int i = 0; i < particles.size(); ++i) {
		QVector3D new_velocity = particles[i].velocity + particles[i].force / particles[i].density * deltaT;
		particles[i].position += (particles[i].velocity + new_velocity) * 0.5 * deltaT;

		// check z coordinates
		if (particles[i].position.z() < 0) {
			particles[i].position.setZ(0.0f);
			new_velocity.setZ(-new_velocity.z());
		} else if (particles[i].position.z() > height_of_container) {
			particles[i].position.setZ(height_of_container);
			new_velocity.setZ(-new_velocity.z());
		}

		// check r
		float r = sqrtf(particles[i].position.x() * particles[i].position.x() + particles[i].position.y() * particles[i].position.y());
		if (r > radius_of_container) {
			particles[i].position.setX(particles[i].position.x() / r * radius_of_container);
			particles[i].position.setY(particles[i].position.y() / r * radius_of_container);

			/*QVector2D v(new_velocity.x(), new_velocity.y());
			QVector2D n(particles[i].position.x(), particles[i].position.y());
			n.normalize();

			v = v * 2 - n * QVector2D::dotProduct(v, n) * 2;
			new_velocity.setX(v.x());
			new_velocity.setY(v.y());
			*/
			new_velocity.setX(0.0f);
			new_velocity.setY(0.0f);
		}

		particles[i].velocity = new_velocity;
	}
	printf("%lf\n", particles[0].position.z());
}

float SPH::W_poly6(float r, float h) {
	return 315.0f / 64 / M_PI / powf(h, 9) * powf(h * h - r * r, 3);
}

QVector3D SPH::dW_spiky(const QVector3D& r, float h) {
	return -45 * (h - r.length()) * (h - r.length()) / M_PI / powf(h, 6) / r.length() * r;
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