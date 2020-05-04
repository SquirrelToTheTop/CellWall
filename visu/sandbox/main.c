void reshape (int w, int h) 
{
	glViewport (0, 0, (GLsizei)w, (GLsizei)h);  
	glMatrixMode (GL_PROJECTION);  
	glLoadIdentity ();  
	gluPerspective (60, (GLfloat)w / (GLfloat)h, 1.0, 500.0);
	glMatrixMode (GL_MODELVIEW);  
}

void renderScene(void)
{
	glClear (GL_COLOR_BUFFER_BIT);  
	glLoadIdentity();
	glRotatef(ry, 0.0, 1.0, 0.0); //rotate about the z axis
	glRotatef(rx-180, 1.0, 0.0, 0.0); //rotate about the y axis  

	float x,y,z;  
	glPointSize(1.0);   
	glBegin(GL_POINTS);  
	for (int i=0;i<480;i++)
	{   
		for (int j=0;j<640;j++)
		{  
			glColor3f(texture[i][j][0]/255, texture[i][j][1]/255, texture[i][j][2]/255);
			x=imgdata[i][j][0];
			y=imgdata[i][j][1];   
			z=imgdata[i][j][2];   
			glVertex3f(x,y,z);   
		}  
	}  
	glEnd();  
	glFlush();  
}

int main(int argc, char **argv)
{
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);  
	glutInitWindowPosition(10,10);
	glutInitWindowSize(640, 480);  
	glutCreateWindow("3D disparity image"); 

	glutDisplayFunc(renderScene);
	glutReshapeFunc(reshape);
	glutIdleFunc(renderScene);
	glutMainLoop();
}
