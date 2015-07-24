#ifndef _EXAMPLE_COMMON_H_
#define _EXAMPLE_COMMON_H_

#ifdef __cplusplus
extern "C" {
#endif
  void drawCircle(float x, float y, float radius, const float* rgba, float angle, char fill, float k_segments/*=16.0f*/);
  void drawPolygon(char fill, float* poly, int count, float x, float y, float* color);
#ifdef __cplusplus
};
#endif

#ifdef EXCOMMON_IMPLEMENTATION
void drawCircle(float x, float y, float radius, const float* rgba, float angle, char fill, float k_segments/*=16.0f*/)
{
  const float k_increment = 2.0f * (float)3.14159f / k_segments;
  float sinInc = sinf(k_increment);
  float cosInc = cosf(k_increment);
  float r1[2]={1.0f, 0.0f};
  float v1[2]={ radius*r1[0], radius*r1[1]};
  float r2[2];
  float v2[2];
  int i;

  glPushMatrix();
  glTranslatef(x,y,.0f);
  glRotatef(angle,0.0f,0.0f,1.0f);
  glColor4fv( (GLfloat*)rgba );
  glBegin(fill?GL_TRIANGLE_FAN:GL_LINE_LOOP);
  for (i = 0; i < k_segments; ++i)
  {
    // Perform rotation to avoid additional trigonometry.
    r2[0] = cosInc * r1[0] - sinInc * r1[1];
    r2[1] = sinInc * r1[0] + cosInc * r1[1];
    v2[0] = radius*r2[0];
    v2[1] = radius*r2[1];

    glVertex2fv( (GLfloat*)v1 );
    glVertex2fv( (GLfloat*)v2 );
    r1[0] = r2[0]; r1[1] = r2[1];
    v1[0] = v2[0]; v1[1] = v2[1];
  }
  glEnd();
  glPopMatrix();
}

void drawPolygon(char fill, float* poly, int count, float x, float y, float* color)
{
  int i;

  glLineWidth(2.0f);
  glPointSize(7.0f);

  glPushMatrix();
  glTranslatef(x,y,0.0f);
  glColor4fv(color);
  if ( fill )
  {
    glBegin(GL_TRIANGLE_FAN);
    for ( i = 0; i < count; ++i )
      glVertex2f(poly[i*2], poly[i*2+1]);
    glEnd();
  }
  else
  {
    glColor4ub(255,255,255,255);
    glBegin(GL_LINE_LOOP);
    for ( i = 0; i < count; ++i )
      glVertex2f(poly[i*2], poly[i*2+1]);
    glEnd();

    glColor4ub(255,255,0,255);
    glBegin(GL_POINTS);
    for ( i = 0; i < count; ++i )
      glVertex2f(poly[i*2], poly[i*2+1]);
    glEnd();
  }
  glPopMatrix();
}

#endif

#endif