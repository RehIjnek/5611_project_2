//Create Window
Camera camera;
PImage img;

String windowTitle = "Swinging Rope";
void setup() {
  size(1500, 900, P3D);
  surface.setTitle(windowTitle);
  camera = new Camera();
  img = loadImage("towel2.jpg");
  initScene();
}
//Obstacle Parameters
Vec3 obstaclePos = new Vec3(200, 170, 100); 
float obstacleRadius = 80;
//Simulation Parameters
float floor = 750;
Vec3 gravity = new Vec3(0,400,0);
float radius = 5;
Vec3 stringTop = new Vec3(200,50,0);
float restLen = 10;
float mass = 1.0; //TRY-IT: How does changing mass affect resting length of the rope?
float k = 200; //TRY-IT: How does changing k affect resting length of the rope?
float kv = 30; //TRY-IT: How big can you make kv?
float dragC = 0.01;
float fluidDens = 0.01;
Vec3 airVel = new Vec3(0, 0, 1000);

//Initial positions and velocities of masses
static int maxNodes = 1000;
Vec3 pos[] = new Vec3[maxNodes];
Vec3 vel[] = new Vec3[maxNodes];
Vec3 acc[] = new Vec3[maxNodes];
boolean torn[] = new boolean[maxNodes];
int numHoriz = 20;
int numVert = 14;
int numNodes = numVert* numHoriz;

float kFric = 30.0;
boolean tearing = false;
void initScene(){
  for (int i = 0; i < numHoriz; i++) {
    for (int j = 0; j < numVert; j++) {
      pos[numVert*i + j] = new Vec3(0,0,0);
      pos[numVert*i + j].x = 100 + restLen*i;
      pos[numVert*i + j].z = restLen*j; //Make each node a little lower
      vel[numVert*i + j] = new Vec3(0,0,0);
      torn[numVert*i + j] = false;
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
  for (int i = 0; i < numHoriz-1; i++){
    for (int j = 0; j < numVert; j++) {
      Vec3 diff = pos[numVert*(i + 1) + j].minus(pos[numVert*i + j]);
      float stringF = -k*(diff.length() - restLen);
    
      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[numVert*i + j], stringDir);
      float projVtop = dot(vel[numVert*(i + 1) + j], stringDir);
      float dampF = -kv*(projVtop - projVbot);
      
      Vec3 force = stringDir.times(stringF+dampF);
      acc[numVert*i + j].add(force.times(-1.0/mass));
      acc[numVert*(i + 1) + j].add(force.times(1.0/mass));  
    }
  }
  
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numVert-1; j++) {
      Vec3 diff = pos[numVert*i + j + 1].minus(pos[numVert*i + j]);
      float stringF = -k*(diff.length() - restLen);
      
      Vec3 stringDir = diff.normalized();
      float projVbot = dot(vel[numVert*i + j], stringDir);
      float projVtop = dot(vel[numVert*i + j + 1], stringDir);
      float dampF = -kv*(projVtop - projVbot);
    
      Vec3 force = stringDir.times(stringF+dampF);
      acc[numVert*i + j].add(force.times(-1.0/mass));
      acc[numVert*i + j + 1].add(force.times(1.0/mass));  
    }
  } 
  
  for (int i = 0; i < numHoriz-1; i++){
    for (int j = 0; j < numVert-1; j++) {
      Vec3 v1 = vel[numVert*i + j];
      Vec3 v2 = vel[numVert*(i + 1) + j];
      Vec3 v3 = vel[numVert*i + j + 1];
      Vec3 v4 = vel[numVert*(i + 1) + j + 1];
      
      Vec3 relVel1 = v1.plus(v3.plus(v2)).times(1.0/3.0).minus(airVel);
      Vec3 relVel2 = v3.plus(v2.plus(v4)).times(1.0/3.0).minus(airVel);
      Vec3 r1 = pos[numVert*i + j];
      Vec3 r2 = pos[numVert*(i + 1) + j];
      Vec3 r3 = pos[numVert*i + j + 1];
      Vec3 r4 = pos[numVert*(i + 1) + j + 1];
      Vec3 nStar1 = cross(r2.minus(r1), r3.minus(r1));
      Vec3 nStar2 = cross(r3.minus(r4), r2.minus(r4));
      Vec3 v2an1 = nStar1.times(relVel1.length() * dot(relVel1, nStar1) / (2 * nStar1.length()));
      Vec3 v2an2 = nStar2.times(relVel2.length() * dot(relVel1, nStar2) / (2 * nStar2.length()));
      Vec3 force1 = v2an1.times(-1.0 / 6.0 * fluidDens * dragC);
      Vec3 force2 = v2an2.times(-1.0 / 6.0 * fluidDens * dragC);
      
      acc[numVert*i + j].add(force1.times(1.0/mass));
      acc[numVert*(i + 1) + j].add(force1.times(1.0/mass));
      acc[numVert*i + j + 1].add(force1.times(1.0/mass));
      acc[numVert*(i + 1) + j].add(force2.times(1.0/mass));
      acc[numVert*i + j + 1].add(force2.times(1.0/mass)); 
      acc[numVert*(i + 1) + j + 1].add(force2.times(1.0/mass));
    }
  } 
  //for (int i = 0; i < numHoriz - 1; i++){
  //  for (int j = 0; j < numVert - 1; j++) {
  //    Vec3 diff = pos[numVert*i + j].minus(pos[numVert*(i + 1) + j + 1]);
  //    float stringF = -k*(diff.length() - restLen); 
    
  //    Vec3 stringDir = diff.normalized();
  //    float projVbot = dot(vel[numVert*i + j], stringDir);
  //    float projVtop = dot(vel[numVert*(i + 1) + j + 1], stringDir);
  //    float dampF = -kv*(projVtop - projVbot);
    
    
  //    Vec3 force = stringDir.times(stringF+dampF);
  //    acc[numVert*i + j].add(force.times(-1.0/mass));
  //    acc[numVert*(i + 1) + j + 1].add(force.times(1.0/mass));  
  //  }
  //} 
  
  //for (int i = 0; i < numHoriz - 1; i++){
  //  for (int j = 0; j < numVert - 1; j++) {
  //    Vec3 diff = pos[numVert*i + j + 1].minus(pos[numVert*(i + 1) + j]);
  //    float stringF = -k*(diff.length() - restLen); 
    
  //    Vec3 stringDir = diff.normalized();
  //    float projVbot = dot(vel[numVert*i + j + 1], stringDir);
  //    float projVtop = dot(vel[numVert*(i + 1) + j], stringDir);
  //    float dampF = -kv*(projVtop - projVbot);
    
  //    Vec3 force = stringDir.times(stringF+dampF);
  //    acc[numVert*i + j + 1].add(force.times(-1.0/mass));
  //    acc[numVert*(i + 1) + j].add(force.times(1.0/mass));  
  //  }
  //} 
  
  // Tearing
  
  for (int i = 0; i < numHoriz; i++){
    for (int j = 1; j < numVert; j++) {
      if(mass * acc[numVert*i + j].length() > 5000) {
        torn[numVert*i + j] = true;
      } 
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
      if (d < obstacleRadius + 0.09) {
        Vec3 n = obstaclePos.minus(pos[numVert*i + j]).times(-1);
        n.normalize();
        Vec3 bounce = n.times(dot(vel[numVert*i + j], n));
        vel[numVert*i + j] = vel[numVert*i + j].minus(bounce.times(1.5));
        pos[numVert*i + j] = pos[numVert*i + j].plus(n.times(1.5 + obstacleRadius - d));
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
  //for (int i = 0; i < numHoriz; i++){
  //  for (int j = 0; j < numVert - 1; j++) {
  //    pushMatrix();
  //    line(pos[numVert*i + j].x, pos[numVert*i + j].y, pos[numVert*i + j].z, pos[numVert*i + j + 1].x, pos[numVert*i + j + 1].y, pos[numVert*i + j + 1].z);
  //    translate(pos[numVert*i + j].x, pos[numVert*i + j].y, pos[numVert*i + j].z);
  //    sphere(radius);
  //    popMatrix();
  //    if(j == numVert - 2) { //Draw the last node of each vertical rope
  //      pushMatrix();
  //      translate(pos[numVert*i + j + 1].x, pos[numVert*i + j + 1].y, pos[numVert*i + j + 1].z);
  //      sphere(radius);
  //      popMatrix();
  //    }
  //  }
  //}
  
  // Draw Ropes Horizontally
  //for (int i = 0; i < numVert; i++){
  //  for (int j = 0; j < numHoriz - 1; j++) {
  //    pushMatrix();
  //    line(pos[numVert*j + i].x, pos[numVert*j + i].y, pos[numVert*j + i].z, pos[numVert*(j + 1) + i].x, pos[numVert*(j + 1) + i].y, pos[numVert*(j + 1) + i].z);
  //    translate(pos[numVert*j + i].x, pos[numVert*j + i].y, pos[numVert*j + i].z);
  //    popMatrix();
  //  }
  //}
  
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
  
  if (paused)
    surface.setTitle(windowTitle + " [PAUSED]");
  else
    surface.setTitle(windowTitle + " "+ nf(frameRate,0,2) + "FPS");
    
  ////create the towel texture
  for (int i = 0; i < numHoriz-1; i++) {
    for (int j = 0; j < numVert-1; j++) { 
      beginShape();
      texture(img);
      normal(1, 0, 1);
      if (!torn[numVert*i + j]) {
        vertex(pos[numVert*i + j].x, pos[numVert*i + j].y, pos[numVert*i + j].z, 0, 0);
      }
      //normal(1, 0, 1);
      if (!torn[numVert*(i + 1) + j]) {
        vertex(pos[numVert*(i + 1) + j].x, pos[numVert*(i + 1) + j].y, pos[numVert*(i + 1) + j].z, img.width, 0);
      }
      //normal(1, 0, 1);
      if (!torn[numVert*(i + 1) + j + 1]) {
        vertex(pos[numVert*(i + 1) + j + 1].x, pos[numVert*(i + 1) + j + 1].y, pos[numVert*(i + 1) + j + 1].z, img.width, img.height);
      }
      //normal(1, 0, 1);
      if (!torn[numVert*i + j + 1]) {
        vertex(pos[numVert*i + j + 1].x, pos[numVert*i + j + 1].y, pos[numVert*i + j + 1].z, 0, img.height);
      }
      endShape();  
    }
  }
  
  stroke(0,0,0);
}

void keyPressed(){
  if (key == ' ') {
    paused = !paused;
  } else if (key == 'r') {
    for (int i = 0; i < numHoriz; i++) {
      for (int j = 0; j < numVert; j++) {
        pos[numVert*i + j] = new Vec3(0,0,0);
        pos[numVert*i + j].x = 100 + restLen*i;
        pos[numVert*i + j].z = restLen*j; //Make each node a little lower
        vel[numVert*i + j] = new Vec3(0,0,0);
        torn[numVert*i + j] = false;
      }
    }
  }
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

Vec3 cross(Vec3 a, Vec3 b) {
    return new Vec3(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x);
}
float dot(Vec3 a, Vec3 b){
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

Vec3 projAB(Vec3 a, Vec3 b){
  return b.times(a.x*b.x + a.y*b.y + a.z*b.z);
}
