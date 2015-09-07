#ifdef __APPLE_CC__
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif

#include <iomanip>
#include "Mesh.h"

int gridX = 600;
int gridY = 600;
int gridZ = 600;

const double fovy = 50.;
const double clipNear = .01;
const double clipFar = 1000.;
double x = 0;
double y = 0;
double z = -2.5;

std::string path = "/Users/rohansawhney/Desktop/developer/C++/geodesics/bunny.obj";

Mesh mesh;
bool success = true;
int vIdx = 0;

void printInstructions()
{
    std::cerr << "' ': draw geodesics\n"
              << "↑/↓: move in/out\n"
              << "w/s: move up/down\n"
              << "a/d: move left/right\n"
              << "escape: exit program\n"
              << std::endl;
}

void init()
{
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glEnable(GL_DEPTH_TEST);
}

void setPhi(double& phi)
{
    for (double d = 0.2; d < 1.2; d += 0.2) {
        if (phi < d) {
            phi = d;
            break;
        }
    }
}

void draw()
{
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (EdgeCIter e = mesh.edges.begin(); e != mesh.edges.end(); e++) {
        double phi = (e->he->vertex->phi + e->he->flip->vertex->phi) / 2.0;
        setPhi(phi);
        
        glColor4f(phi, 0, 0, 0.6);
        glVertex3d(e->he->vertex->position.x(),
                   e->he->vertex->position.y(),
                   e->he->vertex->position.z());
        glVertex3d(e->he->flip->vertex->position.x(),
                   e->he->flip->vertex->position.y(),
                   e->he->flip->vertex->position.z());
        
    }
    glEnd();
    
    glColor4f(1, 1, 1, 0.6);
    glLineWidth(3.0);
    glBegin(GL_LINES);
    HalfEdgeCIter he = mesh.vertices[vIdx].he;
    do {
        glVertex3d(he->vertex->position.x(),
                   he->vertex->position.y(),
                   he->vertex->position.z());
        glVertex3d(he->flip->vertex->position.x(),
                   he->flip->vertex->position.y(),
                   he->flip->vertex->position.z());
        
        he = he->flip->next;
    } while (he != mesh.vertices[vIdx].he);
    glEnd();
}

void display()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    double aspect = (double)viewport[2] / (double)viewport[3];
    gluPerspective(fovy, aspect, clipNear, clipFar);
    
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    
    gluLookAt(0, 1, z, x, y, 0, 0, 1, 0);
    
    if (success) {
        draw();
    }

    glutSwapBuffers();
}

void keyboard(unsigned char key, int x0, int y0)
{
    switch (key) {
        case 27 :
            exit(0);
        case ' ':
            vIdx = rand() % mesh.vertices.size();
            mesh.computeGeodesics(vIdx);
            break;
        case 'a':
            x -= 0.03;
            break;
        case 'd':
            x += 0.03;
            break;
        case 'w':
            y += 0.03;
            break;
        case 's':
            y -= 0.03;
            break;
    }
    
    glutPostRedisplay();
}

void special(int i, int x0, int y0)
{
    switch (i) {
        case GLUT_KEY_UP:
            z += 0.03;
            break;
        case GLUT_KEY_DOWN:
            z -= 0.03;
            break;
        case GLUT_KEY_LEFT:
            break;
        case GLUT_KEY_RIGHT:
            break;
    }
    
    glutPostRedisplay();
}

int main(int argc, char** argv) {
    
    success = mesh.read(path);
    if (success) {
        vIdx = rand() % mesh.vertices.size();
        mesh.computeGeodesics(vIdx);
    }
    
    printInstructions();
    glutInitWindowSize(gridX, gridY);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInit(&argc, argv);
    glutCreateWindow("Geodesics");
    init();
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(special);
    glutMainLoop();
    
    return 0;
}
