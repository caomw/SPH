#include <iostream>
#include "GLWidget3D.h"
#include "MainWindow.h"
#include <GL/GLU.h>
#include "SPH.h"

#define SQR(x)	((x) * (x))

GLWidget3D::GLWidget3D() {

	// set up the camera
	camera.setLookAt(0.0f, 0.0f, 0.0f);
	camera.setYRotation(0);
	camera.setTranslation(0.0f, 0.0f, 2.0f);

	container_width = 2.0;
	container_depth = 0.5;
	container_height = 1.0;
	sph = new SPH(
		0.0457,	// h (radius of smoothing)
		0.02,	// m (mass of particle)
		3,		// k (pressure factor)
		3.5,	// myu (viscosity factor)
		998.29,	// rho_0 (rest density)
		0.00728, // sigma (surface tension coefficient)
		7.065,	// l (surface tension threshold)
		0.6,	// dumping factor c_R
		container_width, container_depth, container_height,
		0.01);	// time step
}

/**
 * This event handler is called when the mouse press events occur.
 */
void GLWidget3D::mousePressEvent(QMouseEvent *e)
{
	lastPos = e->pos();

}

/**
 * This event handler is called when the mouse release events occur.
 */
void GLWidget3D::mouseReleaseEvent(QMouseEvent *e)
{
	updateGL();
}

/**
 * This event handler is called when the mouse move events occur.
 */
void GLWidget3D::mouseMoveEvent(QMouseEvent *e)
{
	float dx = (float)(e->x() - lastPos.x());
	float dy = (float)(e->y() - lastPos.y());
	lastPos = e->pos();

	if (e->buttons() & Qt::LeftButton) {
		camera.changeXRotation(dy);
		camera.changeYRotation(dx);
	} else if (e->buttons() & Qt::RightButton) {
		camera.changeXYZTranslation(0, 0, -dy * camera.dz * 0.02f);
		if (camera.dz < -9000) camera.dz = -9000;
		if (camera.dz > 9000) camera.dz = 9000;
	} else if (e->buttons() & Qt::MidButton) {
		camera.changeXYZTranslation(-dx, dy, 0);
	}

	updateGL();
}

/**
 * This function is called once before the first call to paintGL() or resizeGL().
 */
void GLWidget3D::initializeGL()
{
	glClearColor(0.443, 0.439, 0.458, 0.0);

	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);

	static GLfloat lightPosition[4] = {0.0f, 0.0f, 100.0f, 0.0f};
	glLightfv(GL_LIGHT0, GL_POSITION, lightPosition);

	timer.start(20, this);
}

/**
 * This function is called whenever the widget has been resized.
 */
void GLWidget3D::resizeGL(int width, int height)
{
	height = height?height:1;

	glViewport( 0, 0, (GLint)width, (GLint)height );
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(60, (GLfloat)width/(GLfloat)height, 0.01f, 10000);
	glMatrixMode(GL_MODELVIEW);
}

/**
 * This function is called whenever the widget needs to be painted.
 */
void GLWidget3D::paintGL()
{
	glMatrixMode(GL_MODELVIEW);

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	camera.applyCamTransform();

	drawScene();		
}

/**
 * Draw the scene.
 */
void GLWidget3D::drawScene() {
	sph->draw();

	glLineWidth(1);
	glBegin(GL_LINE_LOOP);
	glColor3f(0, 0, 0);
	glVertex3f(-container_width * 0.5, -container_depth * 0.5, 0);
	glVertex3f(container_width * 0.5, -container_depth * 0.5, 0);
	glVertex3f(container_width * 0.5, container_depth * 0.5, 0);
	glVertex3f(-container_width * 0.5, container_depth * 0.5, 0);
	glEnd();
	glBegin(GL_LINES);
	glVertex3f(-container_width * 0.5, -container_depth * 0.5, 0);
	glVertex3f(-container_width * 0.5, -container_depth * 0.5, container_height);
	glVertex3f(container_width * 0.5, -container_depth * 0.5, 0);
	glVertex3f(container_width * 0.5, -container_depth * 0.5, container_height);
	glVertex3f(container_width * 0.5, container_depth * 0.5, 0);
	glVertex3f(container_width * 0.5, container_depth * 0.5, container_height);
	glVertex3f(-container_width * 0.5, container_depth * 0.5, 0);
	glVertex3f(-container_width * 0.5, container_depth * 0.5, container_height);
	glEnd();
}

QVector2D GLWidget3D::mouseTo2D(int x,int y) {
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	// retrieve the matrices
	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);

	float z;
	glReadPixels(x, (float)viewport[3] - y, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &z);
	
	// unproject the image plane coordinate to the model space
	GLdouble posX, posY, posZ;
	gluUnProject(x, (float)viewport[3] - y, z, modelview, projection, viewport, &posX, &posY, &posZ);

	return QVector2D(posX, posY);
}

void GLWidget3D::timerEvent(QTimerEvent *e) {
    Q_UNUSED(e);

	sph->update();

    updateGL();
}