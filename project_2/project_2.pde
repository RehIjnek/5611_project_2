//Create Window
String windowTitle = "Swinging Rope";
void setup() {
  size(400, 500, P3D);
  surface.setTitle(windowTitle);
  initScene();
}

//Simulation Parameters
float floor = 500;
Vec2 gravity = new Vec2(0,400);
float radius = 5;
Vec2 stringTop = new Vec2(200,50);
float restLen = 10;
float mass = 1.0; //TRY-IT: How does changing mass affect resting length of the rope?
float k = 200; //TRY-IT: How does changing k affect resting length of the rope?
float kv = 30; //TRY-IT: How big can you make kv?

//Initial positions and velocities of masses
static int maxNodes = 100;
Vec2 pos[] = new Vec2[maxNodes];
Vec2 vel[] = new Vec2[maxNodes];
Vec2 acc[] = new Vec2[maxNodes];

int numNodes = 10;
int numHoriz = 5;

float kFric = 30.0;

void initScene(){
  //for (int i = 0; i < numNodes; i++){
  //  pos[i] = new Vec2(0,0);
  //  pos[i].x = stringTop.x;
  //  pos[i].y = stringTop.y + 8*i; //Make each node a little lower
  //  vel[i] = new Vec2(0,0);
  //}
  for (int i = 0; i < numHoriz; i++) {
    for (int j = 0; j < numNodes; j++) {
      pos[10*i + j] = new Vec2(0,0);
      pos[10*i + j].x = 100 + 50*i;
      pos[10*i + j].y = stringTop.y + 8*j; //Make each node a little lower
      vel[10*i + j] = new Vec2(0,0);
    }
  }
}

void update(float dt){

  //Reset accelerations each timestep (momenum only applies to velocity)
  //for (int i = 0; i < numNodes; i++){
  //  acc[i] = new Vec2(0,0);
  //  acc[i].add(gravity);
  //}
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numNodes; j++) {
      acc[10*i + j] = new Vec2(0,0);
      acc[10*i + j].add(gravity);
    }
  }
  
  //Compute (damped) Hooke's law for each spring
  //for (int i = 0; i < numNodes-1; i++){
  //  Vec2 diff = pos[i+1].minus(pos[i]);
  //  float stringF = -k*(diff.length() - restLen);
  //  //println(stringF,diff.length(),restLen);
    
  //  Vec2 stringDir = diff.normalized();
  //  float projVbot = dot(vel[i], stringDir);
  //  float projVtop = dot(vel[i+1], stringDir);
  //  float dampF = -kv*(projVtop - projVbot);
    
  //  float fricF = -kFric*(projVtop - projVbot);
    
  //  Vec2 force = stringDir.times(stringF+dampF+fricF);
  //  acc[i].add(force.times(-1.0/mass));
  //  acc[i+1].add(force.times(1.0/mass));  
  //}
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numNodes-1; j++) {
      Vec2 diff = pos[10*i + j + 1].minus(pos[10*i + j]);
      float stringF = -k*(diff.length() - restLen);
      //println(stringF,diff.length(),restLen);
    
      Vec2 stringDir = diff.normalized();
      float projVbot = dot(vel[10*i + j], stringDir);
      float projVtop = dot(vel[10*i + j + 1], stringDir);
      float dampF = -kv*(projVtop - projVbot);
    
      float fricF = -kFric*(projVtop - projVbot);
    
      Vec2 force = stringDir.times(stringF+dampF+fricF);
      acc[10*i + j].add(force.times(-1.0/mass));
      acc[10*i + j + 1].add(force.times(1.0/mass));  
    }
  }

  //Eulerian integration
  //for (int i = 1; i < numNodes; i++){
  //  vel[i].add(acc[i].times(dt));
  //  pos[i].add(vel[i].times(dt));
  //}
  for (int i = 0; i < numHoriz; i++){
    for (int j = 1; j < numNodes; j++) {
      vel[10*i + j].add(acc[10*i + j].times(dt));
      pos[10*i + j].add(vel[10*i + j].times(dt));
    }
  }
  
  //Collision detection and response
  //for (int i = 0; i < numNodes; i++){
  //  if (pos[i].y+radius > floor){
  //    vel[i].y *= -.9;
  //    pos[i].y = floor - radius;
  //  }
  //}
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numNodes ; j++) {
      if (pos[10*i + j].y+radius > floor){
        vel[10*i + j].y *= -.9;
        pos[10*i + j].y = floor - radius;
      }
    }
  }
  
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
  
  //for (int i = 0; i < numNodes-1; i++){
  //  pushMatrix();
  //  line(pos[i].x,pos[i].y,pos[i+1].x,pos[i+1].y);
  //  translate(pos[i+1].x,pos[i+1].y);
  //  sphere(radius);
  //  popMatrix();
  //}
  for (int i = 0; i < numHoriz; i++){
    for (int j = 0; j < numNodes-1 ; j++) {
      pushMatrix();
      line(pos[10*i + j].x,pos[10*i + j].y,pos[10*i + j + 1].x,pos[10*i + j + 1].y);
      translate(pos[10*i + j + 1].x,pos[10*i + j + 1].y);
      sphere(radius);
      popMatrix();
    }
  }
  
  if (paused)
    surface.setTitle(windowTitle + " [PAUSED]");
  else
    surface.setTitle(windowTitle + " "+ nf(frameRate,0,2) + "FPS");
}

void keyPressed(){
  if (key == ' ')
    paused = !paused;
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
