//Create Window
Camera camera;

String windowTitle = "Swinging Rope";
void setup() {
  size(1500, 900, P3D);
  surface.setTitle(windowTitle);
  camera = new Camera();
  initScene();
}
//Obstacle Parameters
Vec3 obstaclePos = new Vec3(200, 200, 100); 
float obstacleRadius = 80;
//Simulation Parameters
float floor = 500;
Vec3 gravity = new Vec3(0,400,0);
float radius = 5;
Vec3 stringTop = new Vec3(200,50,0);
float restLen = 10;
float mass = 1.0; //TRY-IT: How does changing mass affect resting length of the rope?
float k = 200; //TRY-IT: How does changing k affect resting length of the rope?
float kv = 30; //TRY-IT: How big can you make kv?

//Initial positions and velocities of masses
static int maxNodes = 100;
Vec3 pos[] = new Vec3[maxNodes];
Vec3 vel[] = new Vec3[maxNodes];
Vec3 acc[] = new Vec3[maxNodes];


int numHoriz = 5;
int numVert = 15;
int numNodes = numVert* numHoriz;

float kFric = 30.0;

void initScene(){
  for (int i = 0; i < numHoriz; i++) {
    for (int j = 0; j < numVert; j++) {
      pos[numVert*i + j] = new Vec3(0,0,0);
      pos[numVert*i + j].x = 100 + 50*i;
      pos[numVert*i + j].z = 8*j; //Make each node a little lower
      vel[numVert*i + j] = new Vec3(0,0,0);
    }
  }
}

void update(float dt){

  //Reset accelerations each timestep (momenum only applies to velocity)
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numVert; j++) {
      acc[numVert*i + j] = new Vec3(0,0,0);
      acc[numVert*i + j].add(gravity);
    }
  }
  
  //Compute (damped) Hooke's law for each spring
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numVert-1; j++) {
      Vec3 diff = pos[numVert*i + j + 1].minus(pos[numVert*i + j]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);
    
      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[numVert*i + j], stringDir);
      float projVtop = dot(vel[numVert*i + j + 1], stringDir);
      float dampF = -kv*(projVtop - projVbot);
    
      float fricF = -kFric*(projVtop - projVbot);
    
      Vec3 force = stringDir.times(stringF+dampF+fricF);
      acc[numVert*i + j].add(force.times(-1.0/mass));
      acc[numVert*i + j + 1].add(force.times(1.0/mass));  
    }
  } 
  
  for (int i = 0; i < numHoriz - 1; i++){
    for (int j = 0; j < numVert; j++) {
      Vec3 diff = pos[numVert*(i + 1) + j].minus(pos[numVert*i + j]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);
    
      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[numVert*i + j], stringDir);
      float projVtop = dot(vel[numVert*(i + 1) + j], stringDir);
      float dampF = -kv*(projVtop - projVbot);
    
      float fricF = -kFric*(projVtop - projVbot);
      
      //Vec3 force = stringDir.times(stringF+dampF+fricF);
      Vec3 force = new Vec3(0, 0, 0);
      acc[numVert*i + j].add(force.times(-1.0/mass));
      acc[numVert*(i + 1) + j].add(force.times(1.0/mass));  
    }
  }

  //Eulerian integration
  for (int i = 0; i < numHoriz; i++){
    for (int j = 1; j < numVert; j++) {
      vel[numVert*i + j].add(acc[numVert*i + j].times(dt));
      pos[numVert*i + j].add(vel[numVert*i + j].times(dt));
    }
  }
  
  //Collision detection and response
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numVert ; j++) {
      if (pos[numVert*i + j].y+radius > floor){
        vel[numVert*i + j].y *= -.9;
        pos[numVert*i + j].y = floor - radius;
      }
      
      //sphere collision
      float d = obstaclePos.distanceTo(pos[numVert*i + j]);
      if (d < obstacleRadius+0.09) {
        Vec3 n = obstaclePos.minus(pos[numVert*i + j]).times(-1);
        n.normalize();
        Vec3 bounce = n.times(dot(vel[numVert*i + j], n));
        vel[numVert*i + j] = vel[numVert*i + j].minus(bounce.times(1.5));
        pos[numVert*i + j] = pos[numVert*i + j].plus(n.times(0.1 + obstacleRadius - d));
      }
    }
  }
  
  //update camera
  camera.Update(1.0/frameRate);
}

//Draw the scene: one sphere per mass, one line connecting each pair
boolean paused = true;
void draw() {
  background(255,255,255);
  if (!paused) {
    for (int i = 0; i < 20; i++) {
      update(1/(20*frameRate));
    }
  }
  fill(0,0,0);
  
  // Draw Ropes Vertically
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numVert - 1; j++) {
      pushMatrix();
      line(pos[numVert*i + j].x, pos[numVert*i + j].y, pos[numVert*i + j].z, pos[numVert*i + j + 1].x, pos[numVert*i + j + 1].y, pos[numVert*i + j + 1].z);
      translate(pos[numVert*i + j].x, pos[numVert*i + j].y, pos[numVert*i + j].z);
      sphere(radius);
      popMatrix();
      if(j == numVert - 2) { //Draw the last node of each vertical rope
        pushMatrix();
        translate(pos[numVert*i + j + 1].x, pos[numVert*i + j + 1].y, pos[numVert*i + j + 1].z);
        sphere(radius);
        popMatrix();
      }
    }
  }
  
  // Draw Ropes Horizontally
  for (int i = 0; i < numVert; i++){
    for (int j = 0; j < numHoriz - 1; j++) {
      pushMatrix();
      line(pos[numVert*j + i].x, pos[numVert*j + i].y, pos[numVert*j + i].z, pos[numVert*(j + 1) + i].x, pos[numVert*(j + 1) + i].y, pos[numVert*(j + 1) + i].z);
      translate(pos[numVert*j + i].x, pos[numVert*j + i].y, pos[numVert*j + i].z);
      popMatrix();
    }
  }
  
  pushMatrix();
  fill(0,200,100); 
  specular(120, 120, 180);  //Setup lights… 
  ambientLight(90,90,90);   //More light…
  lightSpecular(255,255,255); 
  shininess(20);  //More light…
  directionalLight(200, 200, 200, -1, 1, -1); //More light…
  translate(obstaclePos.x, obstaclePos.y, obstaclePos.z);
    noStroke();
  sphere(obstacleRadius);   //Draw sphere

  popMatrix();
  
  stroke(0,0,0);
  if (paused)
    surface.setTitle(windowTitle + " [PAUSED]");
  else
    surface.setTitle(windowTitle + " "+ nf(frameRate,0,2) + "FPS");
}

void keyPressed(){
  if (key == ' ')
    paused = !paused;
  camera.HandleKeyPressed();
}

void keyReleased() {
  camera.HandleKeyReleased();
}

///////////////////
// Vec2D Library
///////////////////

public class Vec2 {
  public float x, y;
  
  public Vec2(float x, float y){
    this.x = x;
    this.y = y;
  }
  
  public String toString(){
    return "(" + x+ ", " + y +")";
  }
  
  public float length(){
    return sqrt(x*x+y*y);
  }
  
  public float lengthSqr(){
    return x*x+y*y;
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
  
  public void normalize(){
    float magnitude = sqrt(x*x + y*y);
    x /= magnitude;
    y /= magnitude;
  }
  
  public Vec2 normalized(){
    float magnitude = sqrt(x*x + y*y);
    return new Vec2(x/magnitude, y/magnitude);
  }
  
  public void clampToLength(float maxL){
    float magnitude = sqrt(x*x + y*y);
    if (magnitude > maxL){
      x *= maxL/magnitude;
      y *= maxL/magnitude;
    }
  }
  
  public void setToLength(float newL){
    float magnitude = sqrt(x*x + y*y);
    x *= newL/magnitude;
    y *= newL/magnitude;
  }
  
  public float distanceTo(Vec2 rhs){
    float dx = rhs.x - x;
    float dy = rhs.y - y;
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

Vec2 projAB(Vec2 a, Vec2 b){
  return b.times(a.x*b.x + a.y*b.y);
}

///////////////////
// Vec3D Library
///////////////////

public class Vec3 {
  public float x, y, z;
  
  public Vec3(float x, float y, float z){
    this.x = x;
    this.y = y;
    this.z = z;
  }
  
  public String toString(){
    return "(" + x + ", " + y + ", " + z + ")";
  }
  
  public float length(){
    return sqrt(x*x+y*y+z*z);
  }
  
  public float lengthSqr(){
    return x*x+y*y+z*z;
  }
  
  public Vec3 plus(Vec3 rhs){
    return new Vec3(x+rhs.x, y+rhs.y, z+rhs.z);
  }
  
  public void add(Vec3 rhs){
    x += rhs.x;
    y += rhs.y;
    z += rhs.z;
  }
  
  public Vec3 minus(Vec3 rhs){
    return new Vec3(x-rhs.x, y-rhs.y, z-rhs.z);
  }
  
  public void subtract(Vec3 rhs){
    x -= rhs.x;
    y -= rhs.y;
    z -= rhs.z;
  }
  
  public Vec3 times(float rhs){
    return new Vec3(x*rhs, y*rhs, z*rhs);
  }
  
  public void mul(float rhs){
    x *= rhs;
    y *= rhs;
    z *= rhs;
  }
  
  public void normalize(){
    float magnitude = sqrt(x*x + y*y + z*z);
    x /= magnitude;
    y /= magnitude;
    z /= magnitude;
  }
  
  public Vec3 normalized(){
    float magnitude = sqrt(x*x + y*y + z*z);
    return new Vec3(x/magnitude, y/magnitude, z/magnitude);
  }
  
  public void clampToLength(float maxL){
    float magnitude = sqrt(x*x + y*y + z*z);
    if (magnitude > maxL){
      x *= maxL/magnitude;
      y *= maxL/magnitude;
      z *= maxL/magnitude;
    }
  }
  
  public void setToLength(float newL){
    float magnitude = sqrt(x*x + y*y + z*z);
    x *= newL/magnitude;
    y *= newL/magnitude;
    z *= newL/magnitude;
  }
  
  public float distanceTo(Vec3 rhs){
    float dx = rhs.x - x;
    float dy = rhs.y - y;
    float dz = rhs.z - z;
    return sqrt(dx*dx + dy*dy + dz*dz);
  }
  
}

float dot(Vec3 a, Vec3 b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec3 projAB(Vec3 a, Vec3 b){
  return b.times(a.x*b.x + a.y*b.y + a.z*b.z);
}
