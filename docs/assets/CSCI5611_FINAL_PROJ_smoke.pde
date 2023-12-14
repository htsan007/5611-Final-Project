    
//Variables 
int N = 250-2; //smaller runs better
//int N = 498; //grid size (num squares) for 500 pixel dimensions (minus boundary square 499-1)
float h = 1/N;
float diff = 0.00001;//0.00001; 
int size = (N + 2) * (N + 2); //array size
float visc = 0; //0.0001;
boolean paused = false;
float u[] = new float[size];
float v[] = new float[size];
float u_prev[] = new float[size];
float v_prev[] = new float[size];
float dens[] = new float[size];
float dens_prev[] = new float[size];

//advect oxy/fuel too? How to consume and turn into heat?

//Let heat be on 1-0 gradient hottest fire to disapating smoke
//decrease T by dt*cooling_rate (based on fire or smoke, set by diff. thresholds)
//advect T along vel


//set thresholds for color, too small heat don't show, small be black, 
//larger gray smoke, red, orange, yellow, white, blue (closer to source is hotter)

//make separate updateFire function? 


int IX(int i, int j){ //for defintion indexing array
    return (i + (N + 2) * j);
}

void SWAP(float x0[], float x[]){
  float temp[] = new float[x0.length];
  arrayCopy(x0,temp);
  arrayCopy(x,x0);
  arrayCopy(temp,x);
  //temp = x0;
  //x0 = x;
  //x = temp;
  
}

void setup(){
  //size(500,500,P2D);
  size(250, 250, P2D); //P2D for 2d
  //frameRate(60);
  surface.setTitle("CSCI5611 Final Project SMOKE Sim");
}

void keyPressed(){
  if (key == ' '){
    paused = !paused;
  }
}

void draw() {
  float dt = 1.0 / 40;
  //for(int i = 0; i < 10; i++){
    update_physics(dt);
  
  for(int i = 0; i < N; i++){
    for(int j = 0; j < N; j++){
     
      float col = dens[IX(i,j)]; //compute smoke color from density
      
      color c = color(col,col,col);
      set(i, j, c);
      
    }
  }
}
void update_physics(float dt){
  int i,j;
  if(!paused){
    for(i = 120; i <= 140; i++){  //smaller scene version
       for(j = 220; j <= 248; j+=3){
        dens_prev[IX(i,j)] = 555;
        dens_prev[IX(i,j+1)] = 555;
        dens_prev[IX(i,j+2)] = 555;
        v_prev[IX(i,j+2)] = -3;
        v_prev[IX(i,j+1)] = -4;
        v_prev[IX(i,j)] = -6;
        
      }
    }
  }
  else {
    for(i = 120; i <= 140; i++){  //smaller scene version
     for(j = 220; j <= 248; j+=3){
       dens_prev[IX(i,j)] = 0;
       dens_prev[IX(i,j+1)] = 0;
       dens_prev[IX(i,j+2)] = 0;
      v_prev[IX(i,j)] = 0;
      v_prev[IX(i,j+1)] = 0;
      v_prev[IX(i,j+2)] = 0;
      }
    }
  }

  
  vel_step(N, u, v, u_prev, v_prev, visc, dt);
  dens_step(N, dens, dens_prev, u, v, diff, dt);
  for(int p = 0; p < size; p++){ //fade smoke
    dens[p] *= 0.975;
  }
 
}

void add_source (int N, float x[], float s[], float dt){
  //size = (N+2)*(N+2);
  int i;
  for(i = 0; i < size; i++){
    x[i] += dt * s[i];
  }
}

void set_bnd (int N, int b, float x[]){
  int i;
  for(i = 1; i<=N; i++){
    if( b==1){
      x[IX(0,i)] = -x[IX(1,i)];
      x[IX(N+1,i)] = -x[IX(N,i)];
    }
    else {
      x[IX(0,i)] = x[IX(1,i)];
      x[IX(N+1,i)] = x[IX(N,i)];
    }
    if( b==2){
      x[IX(i,0)] = -x[IX(i,1)];
      x[IX(i,N+1)] = -x[IX(i,N)];
    }
    else {
      x[IX(i,0)] = x[IX(i,1)];
      x[IX(i,N+1)] = x[IX(i,N)];
    }
    x[IX(0 ,0 )] = 0.5 * (x[IX(1,0 )] + x[IX(0 ,1)]);
    x[IX(0 ,N+1)] = 0.5 * (x[IX(1,N+1)] + x[IX(0 ,N )]);
    x[IX(N+1,0 )] = 0.5 * (x[IX(N,0 )] + x[IX(N+1,1)]);
    x[IX(N+1,N+1)] = 0.5 * (x[IX(N,N+1)] + x[IX(N+1,N )]);
  }
}

void diffuse (int N, int b, float x[], float x0[], float diff, float dt){
  float a = dt * diff * N * N;
  int i,j,k;
  for(k = 0; k < 20; k++){
    for(i = 1; i <= N; i++){
      for(j = 1; j <= N; j++){
        x[IX(i,j)] = (x0[IX(i,j)] + a * (x[IX(i-1, j)] + x[IX(i+1, j)] 
                                      + x[IX(i, j-1)] + x[IX(i, j+1)])) / (1 + 4 * a);
      }
    }
    set_bnd(N ,b ,x);
  }
}

void advect (int N, int b, float d[], float d0[], float u[], float v[], float dt){
  
  int i, j, i0, j0, i1, j1;
  float x, y, s0, t0, s1, t1, dt0;
  
  dt0 = dt*N;
  for( i = 1; i<= N; i++){
    for( j = 1; j<= N; j++){
      x = i - dt0 * u[IX(i,j)];
      y = j - dt0*v[IX(i,j)];
      if(x < 0.5){
        x = 0.5;
      }
      if( x > N + 0.5){ //swap to if-else if??
        x = N + 0.5;
      }
      i0 = (int)x;
      i1 = i0 + 1;
      
      if(y < 0.5){
        y= 0.5;
      }
      if(y > N + 0.5){
        y = N + 0.5;
      }
      j0 = (int)y;
      j1 = j0 + 1;
      
      s1 = x - i0;
      s0 = 1 - s1;
      t1 = y - j0;
      t0 = 1 - t1;
      
      d[IX(i,j)] = s0 * (t0 * d0[IX(i0,j0)] + t1 * d0[IX(i0,j1)]) 
                + s1 *(t0 * d0[IX(i1,j0)] + t1 * d0[IX(i1,j1)]);
    }
  }
  set_bnd(N,b,d);
}

void project( int N, float u[], float v[], float p[], float div[]){
  int i, j, k;
  float h = 1.0/N;
  
  for( i = 1; i <= N; i++){
    for( j = 1; j <= N; j++){
      div[IX(i,j)] = -0.5 * h * (u[IX(i+1,j)] - u[IX(i-1,j)] +
                                  v[IX(i,j+1)] - v[IX(i,j-1)]);
      p[IX(i,j)] = 0;
    }
  }
  set_bnd(N, 0 ,div);
  set_bnd(N, 0 ,p);
  
  for(k = 0; k < 20; k++){
    for(i = 1; i <= N; i++){
      for(j = 1; j <= N; j++){
        p[IX(i,j)] = (div[IX(i,j)] + p[IX(i-1,j)] + p[IX(i+1,j)] +
                               p[IX(i,j-1)] + p[IX(i,j+1)]) / 4;
      }
    }
    set_bnd(N ,0 ,p );
  }
  
   for(i = 1; i <= N; i++){
      for(j = 1; j <= N; j++){
        u[IX(i,j)] -= 0.5 * (p[IX(i+1,j)] - p[IX(i-1,j)]) / h;
        v[IX(i,j)] -= 0.5 * (p[IX(i,j+1)] - p[IX(i,j-1)]) / h;
      }
   }
   set_bnd(N, 1, u);
   set_bnd(N, 2, v);
}

void dens_step( int N, float x[], float x0[], float u[], float v[], float diff, float dt){
  add_source(N, x, x0, dt);
  SWAP(x0,x);
  diffuse(N, 0, x, x0, diff, dt);
  SWAP(x0,x);
  advect(N, 0, x, x0, u, v, dt);
}

void vel_step(int N, float u[], float v[], float u0[], float v0[], float visc, float dt){
  add_source(N, u, u0, dt);
  add_source(N, v, v0, dt);
  SWAP(u0, u);
  diffuse(N, 1, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(N, 2, v, v0, visc, dt);
  
  project(N, u, v, u0, v0);
  SWAP(u0, u);
  SWAP(v0, v);
  advect(N, 1, u, u0, u0, v0, dt);
  advect(N, 2, v, v0, u0, v0, dt);
  project(N, u, v, u0, v0);
  
}











//---------------
//Vec 2 Library
//---------------

//2DVector library
public class Vec2 {
  public float x, y;
  
  public Vec2(float x, float y){
    this.x = x;
    this.y = y;
  }
  
  public String toString(){
    return "(" + x+ "," + y +")";
  }
  
  public float length(){
    if((x==0) && (y==0)){
      return 0;
    }
    else {
      return sqrt(x*x+y*y);
    }
  }
  
  public float lengthSqr(){
    return x*x + y*y;
  }
  
  
  public Vec2 plus(Vec2 rhs){
    return new Vec2(x+rhs.x, y+rhs.y);
  }
  
  public void add(Vec2 rhs){
    x += rhs.x;
    y += rhs.y;
  }
  
  public Vec2 minus(Vec2 rhs){
    return new Vec2(x-rhs.x, y-rhs.y);
  }
  
  public void subtract(Vec2 rhs){
    x -= rhs.x;
    y -= rhs.y;
  }
  
  public Vec2 times(float rhs){
    return new Vec2(x*rhs, y*rhs);
  }
  
  public void mul(float rhs){
    x *= rhs;
    y *= rhs;
  }
  
  public void clampToLength(float maxL){
    if((x==0) && (y==0)){
     return;
    }
    float magnitude = sqrt(x*x + y*y);
    if (magnitude > maxL){
      x *= maxL/magnitude;
      y *= maxL/magnitude;
    }

  }
  
  
  public void setToLength(float newL){
    if((x==0) && (y==0)){
     return;
    }
    float magnitude = sqrt(x*x + y*y);
    x *= newL/magnitude;
    y *= newL/magnitude;
  }
  
  public void normalize(){
    if((x == 0) && (y == 0)){
      return;
    }
    else {
      float magnitude = sqrt(x*x + y*y);
      x /= magnitude;
      y /= magnitude;
    }
  }
  
  public Vec2 normalized(){
    if((x==0) && (y==0)){
     return new Vec2(0,0);
    }
    float magnitude = sqrt(x*x + y*y);
    return new Vec2(x/magnitude, y/magnitude);
  }
  
  public float distanceTo(Vec2 rhs){
    float dx = rhs.x - x;
    float dy = rhs.y - y;
    if((dx == 0) && (dy == 0)){
      return 0;
    }
    return sqrt(dx*dx + dy*dy);
  }
}

Vec2 interpolate(Vec2 a, Vec2 b, float t){
  return a.plus((b.minus(a)).times(t));
}

float interpolate(float a, float b, float t){
  return a + ((b-a)*t);
}

float dot(Vec2 a, Vec2 b){
  return a.x*b.x + a.y*b.y;
}
float vecCross(Vec2 a, Vec2 b){
  return (a.x*b.y - a.y*b.x); //ad - bc
}
Vec2 projAB(Vec2 a, Vec2 b){
  return b.times(a.x*b.x + a.y*b.y);
}

Vec2 perpendicular(Vec2 a) {
  return new Vec2(-a.y, a.x);
}
