#pragma once

#include <QGLWidget>
#include <QMouseEvent>
#include <QKeyEvent>
#include "Camera.h"
#include <QVector3D>
#include <vector>
#include <QBasicTimer>

using namespace std;

class MainWindow;
class SPH;

class GLWidget3D : public QGLWidget {
private:
	QBasicTimer timer;
	Camera camera;
	QPoint lastPos;
	SPH* sph;

public:
	GLWidget3D();
	void drawScene();
	QVector2D mouseTo2D(int x,int y);

protected:
	void initializeGL();
	void resizeGL(int width, int height);
	void paintGL();    
	void mousePressEvent(QMouseEvent *e);
	void mouseMoveEvent(QMouseEvent *e);
	void mouseReleaseEvent(QMouseEvent *e);
	void timerEvent(QTimerEvent *e);
};

