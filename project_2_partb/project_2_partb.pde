
//  3. Add a second box to the scene! Break this into a few steps
//       1) Allow the new box to interact with the 4 walls and the click force
//       2) Detect if the corner of one box is inside the other box
//       3) Use the resources mentioned in resolveCollision() to update that function for object-object collisions
//     If you get this working, you'll have the basics of an impressive 2D physics engine
//  4. Other ideas: Allow other shapes besides boxes, add friction, add gravity

void setup(){
  size(1200,800,P3D);
  shelf = loadImage("shelf.jpeg");
  book = loadImage("book.jpg");
}
PImage shelf;
PImage book;

//Set inital conditions
float w = 50;
float h = 200;
float box_bounce = 1; //Coef. of restitution

float mass = 1;                         //Resistance to change in momentum/velocity
float rot_inertia = mass*(w*w+h*h)/12;  //Resistance to change in angular momentum/angular velocity

Vec2 momentum = new Vec2(0,0);          //Speed the box is translating (derivative of position)
float angular_momentum = 0;             //Speed the box is rotating (derivative of angle)

Vec2 center = new Vec2(400,400);        //Current position of center of mass
float angle = 1.57; /*radians*/          //Current rotation amount (orientation)

Vec2 total_force = new Vec2(0,0);       //Forces change position (center of mass)
float total_torque = 0;                 //Torques change orientation (angle)

Vec2 p1,p2,p3,p4;                       //4 corners of the box -- computed in updateCornerPositions()

//Set inital conditions
Vec2 momentum2 = new Vec2(0,0);          //Speed the box is translating (derivative of position)
float angular_momentum2 = 0;             //Speed the box is rotating (derivative of angle)

Vec2 center2 = new Vec2(800,400);        //Current position of center of mass
float angle2 = 0; /*radians*/          //Current rotation amount (orientation)

Vec2 total_force2 = new Vec2(0,0);       //Forces change position (center of mass)
float total_torque2 = 0;                 //Torques change orientation (angle)

Vec2 q1,q2,q3,q4;                       //4 corners of the box -- computed in updateCornerPositions()

float arrow_angle = 0;

//----------
// Physics Functions
void apply_force(Vec2 force, Vec2 applied_position, int boxNum){
  if (boxNum == 1) {
    total_force.add(force);
    Vec2 displacement = applied_position.minus(center);
    total_torque += cross(displacement, force);
  }
  if (boxNum == 2) {
    total_force2.add(force);
    Vec2 displacement2 = applied_position.minus(center2);
    total_torque2 += cross(displacement2, force);
  }
}

void update_physics(float dt){
  //Update center of mass
  momentum.add(total_force.times(dt));     //Linear Momentum = Force * time
  Vec2 box_vel = momentum.times(1.0/mass); //Velocity = Momentum / mass
  center.add(box_vel.times(dt));           //Position += Vel * time
  
  angular_momentum += total_torque * dt;
  float angular_vel = angular_momentum / rot_inertia;
  angle += angular_vel * dt;
  
  //Reset forces and torques
  total_force = new Vec2(0,0); //Set forces to 0 after they've been applied
  total_torque = 0; //Set torques to 0 after the forces have been applied
  
  //Update center of mass
  momentum2.add(total_force2.times(dt));     //Linear Momentum = Force * time
  Vec2 box_vel2 = momentum2.times(1.0/mass); //Velocity = Momentum / mass
  center2.add(box_vel2.times(dt));           //Position += Vel * time
  
  angular_momentum2 += total_torque2 * dt;
  float angular_vel2 = angular_momentum2 / rot_inertia;
  angle2 += angular_vel2 * dt;
  
  //Reset forces and torques
  total_force2 = new Vec2(0,0); //Set forces to 0 after they've been applied 
  total_torque2 = 0; //Set torques to 0 after the forces have been applied
}


class ColideInfo{
  public boolean hit = false;
  public Vec2 hitPoint = new Vec2(0,0);
  public Vec2 objectNormal =  new Vec2(0,0);
  public int boxNum = 0;
}


void updateCornerPositions(){
  Vec2 right = new Vec2(cos(angle),sin(angle)).times(w/2);
  Vec2 up = new Vec2(-sin(angle),cos(angle)).times(-h/2);
  p1 = center.plus(right).plus(up); //bottom right corner
  p2 = center.plus(right).minus(up); //top right corner
  p3 = center.minus(right).plus(up); //bottom left corner
  p4 = center.minus(right).minus(up); //top left corner
  
  Vec2 right2 = new Vec2(cos(angle2),sin(angle2)).times(w/2);
  Vec2 up2 = new Vec2(-sin(angle2),cos(angle2)).times(-h/2);
  q1 = center2.plus(right2).plus(up2); //bottom right corner
  q2 = center2.plus(right2).minus(up2); //top right corner
  q3 = center2.minus(right2).plus(up2); //bottom left corner
  q4 = center2.minus(right2).minus(up2); //top left corner
}

ColideInfo collisionTest1(){
  updateCornerPositions(); //Compute the 4 corners: p1,p2,p3,p4
  //We only check if the corners collide
  
  ColideInfo info = new ColideInfo();
  //check if it hits right wall
  if (p1.x > width){
    info.hitPoint = p1;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 1;
  }
  if (p2.x > width){
    info.hitPoint = p2;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 1;
  }
  if (p3.x > width){
    info.hitPoint = p3;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 1;
  }
  if (p4.x > width){
    info.hitPoint = p4;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 1;
  }
  //TODO: Test the 4 corners against the left wall
  if (p1.x < 0){
    info.hitPoint = p1;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 1;
  }
  if (p2.x < 0){
    info.hitPoint = p2;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 1;
  }
  if (p3.x < 0){
    info.hitPoint = p3;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 1;
  }
  if (p4.x < 0){
    info.hitPoint = p4;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 1;
  }
  //check if it hits bottom wall
  if (p1.y > height){
    info.hitPoint = p1;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 1;
  }
  if (p2.y > height){
    info.hitPoint = p2;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 1;
  }
  if (p3.y > height){
    info.hitPoint = p3;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 1;
  }
  if (p4.y > height){
    info.hitPoint = p4;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 1;
  }
  //check if it hits top wall
  if (p1.y < 0){
    info.hitPoint = p1;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 1;
  }
  if (p2.y < 0){
    info.hitPoint = p2;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 1;
  }
  if (p3.y < 0){
    info.hitPoint = p3;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 1;
  }
  if (p4.y < 0){
    info.hitPoint = p4;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 1;
  }
  
  return info;
}

ColideInfo collisionTest2(){
  updateCornerPositions(); //Compute the 4 corners: q1,q2,q3,q4
  //We only check if the corners collide
  
  ColideInfo info = new ColideInfo();
  
  //check if it hits right wall
  if (q1.x > width){
    info.hitPoint = q1;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 2;
  }
  if (q2.x > width){
    info.hitPoint = q2;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 2;
  }
  if (q3.x > width){
    info.hitPoint = q3;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 2;
  }
  if (q4.x > width){
    info.hitPoint = q4;
    info.hit = true;
    info.objectNormal = new Vec2(-1,0);
    info.boxNum = 2;
  }
  //TODO: Test the 4 corners against the left wall
  if (q1.x < 0){
    info.hitPoint = q1;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 2;
  }
  if (q2.x < 0){
    info.hitPoint = q2;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 2;
  }
  if (q3.x < 0){
    info.hitPoint = q3;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 2;
  }
  if (q4.x < 0){
    info.hitPoint = q4;
    info.hit = true;
    info.objectNormal = new Vec2(1,0);
    info.boxNum = 2;
  }
  //check if it hits bottom wall
  if (q1.y > height){
    info.hitPoint = q1;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 2;
  }
  if (q2.y > height){
    info.hitPoint = q2;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 2;
  }
  if (q3.y > height){
    info.hitPoint = q3;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 2;
  }
  if (q4.y > height){
    info.hitPoint = q4;
    info.hit = true;
    info.objectNormal = new Vec2(0,-1);
    info.boxNum = 2;
  }
  //check if it hits top wall
  if (q1.y < 0){
    info.hitPoint = q1;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 2;
  }
  if (q2.y < 0){
    info.hitPoint = q2;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 2;
  }
  if (q3.y < 0){
    info.hitPoint = q3;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 2;
  }
  if (q4.y < 0){
    info.hitPoint = q4;
    info.hit = true;
    info.objectNormal = new Vec2(0,1);
    info.boxNum = 2;
  }
  
  return info;
}

//object to object collision test
ColideInfo objectCollisionTest(){
  updateCornerPositions(); //Compute the 4 corners: p1,p2,p3,p4
  //We only check if the corners collide
  
  ColideInfo info = new ColideInfo();
  
  //object to object collision
  if (point_in_box(p1, center2, w , h, angle2)) {
    info.hitPoint = p1;
    info.hit = true; 
    info.boxNum = 1;
  }
  if (point_in_box(p2, center2, w , h, angle2)) {
    info.hitPoint = p2;
    info.hit = true;  
    info.boxNum = 1;
  }
  if (point_in_box(p3, center2, w , h, angle2)) {
    info.hitPoint = p3;
    info.hit = true; 
    info.boxNum = 1;
  }
  if (point_in_box(p4, center2, w , h, angle2)) {
    info.hitPoint = p4;
    info.hit = true;  
    info.boxNum = 1;
  }
  if (point_in_box(q1, center, w , h, angle)) {
    info.hitPoint = q1;
    info.hit = true; 
    info.boxNum = 2;
  }
  if (point_in_box(q2, center, w , h, angle)) {
    info.hitPoint = q2;
    info.hit = true;  
    info.boxNum = 2;
  }
  if (point_in_box(q3, center, w , h, angle)) {
    info.hitPoint = q3;
    info.hit = true; 
    info.boxNum = 2;
  }
  if (point_in_box(q4, center, w , h, angle)) {
    info.hitPoint = q4;
    info.hit = true;
    info.boxNum = 2;
  }
  
  //set object normals
  if (info.boxNum == 1) {
    Vec2 collisionPoint = info.hitPoint.times(1);
    Vec2 relative_pos = collisionPoint.minus(center2);
    Vec2 box_right = new Vec2(cos(angle2),sin(angle2));
    Vec2 box_up = new Vec2(sin(angle2),-cos(angle2));
    Vec2 box_left = box_right.times(-1);
    Vec2 box_down = box_up.times(-1);
    float point_right = dot(relative_pos,box_right);
    float point_up = dot(relative_pos,box_up);
    float distEdges[] = new float[4];
    distEdges[0] = abs(point_right + w/2);
    distEdges[1] = abs(point_right - w/2);
    distEdges[2] = abs(point_up + h/2);
    distEdges[3] = abs(point_up - h/2);
    float min = distEdges[0];
    float minIndex = 0;
    for (int i = 1; i < distEdges.length; i++) { 
      if (distEdges[i] < min) {
        min = distEdges[i];
        minIndex = i;
      }
    }
    if (minIndex == 0) {
      info.objectNormal = box_left;
    }
    if (minIndex == 1) {
      info.objectNormal = box_right;
    }
    if (minIndex == 2) {
      info.objectNormal = box_down;
    }
    if (minIndex == 3) {
      info.objectNormal = box_up;
    }   
  }
  if (info.boxNum == 2) {
    Vec2 collisionPoint = info.hitPoint.times(1);
    Vec2 relative_pos = collisionPoint.minus(center);
    Vec2 box_right = new Vec2(cos(angle),sin(angle));
    Vec2 box_up = new Vec2(sin(angle),-cos(angle));
    Vec2 box_left = box_right.times(-1);
    Vec2 box_down = box_up.times(-1);
    float point_right = dot(relative_pos,box_right);
    float point_up = dot(relative_pos,box_up);
    float distEdges[] = new float[4];
    distEdges[0] = abs(point_right + w/2);
    distEdges[1] = abs(point_right - w/2);
    distEdges[2] = abs(point_up + h/2);
    distEdges[3] = abs(point_up - h/2);
    float min = distEdges[0];
    float minIndex = 0;
    for (int i = 1; i < distEdges.length; i++) { 
        if (distEdges[i] < min) {
            min = distEdges[i];
            minIndex = i;
        }
    }
    if (minIndex == 0) {
      info.objectNormal = box_left;
    }
    if (minIndex == 1) {
      info.objectNormal = box_right;
    }
    if (minIndex == 2) {
      info.objectNormal = box_down;
    }
    if (minIndex == 3) {
      info.objectNormal = box_up;
    }   
  }
  
  return info;
}

//Updates momentum & angular_momentum based on collision using an impulse based method
//This method assumes you hit an immovable obstacle which simplifies the math
// see Eqn 8-18 of here: https://www.cs.cmu.edu/~baraff/sigcourse/notesd2.pdf
// or Eqn 9 here: http://www.chrishecker.com/images/e/e7/Gdmphys3.pdf
//for obstacle-obstacle collisions.
void resolveCollision(Vec2 hit_point, Vec2 hit_normal, float dt, int boxNum){
  if (boxNum == 1) {
    center.add(hit_normal.times(2));
    Vec2 r = hit_point.minus(center);
    Vec2 r_perp = perpendicular(r);
    Vec2 object_vel = momentum.times(1/mass);
    float object_angular_speed = angular_momentum/rot_inertia;
    Vec2 point_vel = object_vel.plus(r_perp.times(object_angular_speed));
    float j = -(1+box_bounce)*dot(point_vel,hit_normal);
    j /= (1/mass + pow(dot(r_perp,hit_normal),2)/rot_inertia);
 
    Vec2 impulse = hit_normal.times(j);
    momentum.add(impulse);
    angular_momentum += dot(r_perp,impulse);
    update_physics(1.01*dt);
  }
  
  if (boxNum == 2) {
    center2.add(hit_normal.times(2));
    Vec2 r2 = hit_point.minus(center2);
    Vec2 r_perp2 = perpendicular(r2);
    Vec2 object_vel2 = momentum2.times(1/mass);
    float object_angular_speed2 = angular_momentum2/rot_inertia;
    Vec2 point_vel2 = object_vel2.plus(r_perp2.times(object_angular_speed2));
    float j = -(1+box_bounce)*dot(point_vel2,hit_normal);
    j /= (1/mass + pow(dot(r_perp2,hit_normal),2)/rot_inertia);
 
    Vec2 impulse2 = hit_normal.times(j);
    momentum2.add(impulse2);
    angular_momentum2 += dot(r_perp2,impulse2);
    update_physics(1.01*dt);
  }
}

void resolveObjectCollision(Vec2 hit_point, Vec2 hit_normal, float dt, int boxNum){
  Vec2 rAP = new Vec2 (0,0);
  Vec2 rBP = new Vec2 (0,0);
  Vec2 vAB = new Vec2 (0,0);
  if (boxNum == 1) {
    center.add(hit_normal.times(2));
    rAP = perpendicular(hit_point.minus(center));    
    rBP = perpendicular(hit_point.minus(center2));
    vAB = momentum.times(1.0/mass).minus(momentum2.times(1.0/mass));
    float j = -(1 + box_bounce) * dot(vAB, hit_normal);
    j /= (dot(hit_normal, hit_normal) * (1.0/mass + 1.0/mass)) + (pow(dot(rAP, hit_normal), 2) / rot_inertia) + (pow(dot(rBP, hit_normal), 2) / rot_inertia);
    
    Vec2 impulse = hit_normal.times(j);
    
    momentum.add(impulse);
    angular_momentum += dot(rAP, impulse);
    
    float j2 = -j;
    Vec2 impulse2 = hit_normal.times(j2);
    
    momentum2.add(impulse2);
    angular_momentum2 += dot(rBP, impulse2);
    
    update_physics(1.01*dt);
  }
  
  if (boxNum == 2) {    
    center2.add(hit_normal.times(2));
    rAP = perpendicular(hit_point.minus(center2));
    rBP = perpendicular(hit_point.minus(center));
    vAB = momentum2.times(1.0/mass).minus(momentum.times(1.0/mass));
    float j = -(1 + box_bounce) * dot(vAB, hit_normal);
    j /= (dot(hit_normal, hit_normal) * (1.0/mass + 1.0/mass)) + (pow(dot(rAP, hit_normal), 2) / rot_inertia) + (pow(dot(rBP, hit_normal), 2) / rot_inertia);
    
    Vec2 impulse = hit_normal.times(j);
    
    momentum2.add(impulse);
    angular_momentum2 += dot(rAP, impulse);
    
    float j2 = -j;
    Vec2 impulse2 = hit_normal.times(j2);
    
    momentum.add(impulse2);
    angular_momentum += dot(rBP, impulse2);
    
    update_physics(1.01*dt);
  }
}

void draw(){
  float dt = 1/(2 * frameRate);
  update_physics(dt);
  
  boolean clicked_box = mousePressed && point_in_box(new Vec2(mouseX, mouseY),center,w,h,angle);
  boolean clicked_box2 = mousePressed && point_in_box(new Vec2(mouseX, mouseY),center2,w,h,angle2);
  
  if (clicked_box) {
    Vec2 force = new Vec2(1,0).times(100);
    if (arrow_angle == 0) {
      force = new Vec2(1,0).times(100);
    }
    if (arrow_angle == 90) {
      force = new Vec2(0,1).times(100);
    }
    if (arrow_angle == 180) {
      force = new Vec2(-1,0).times(100);
    }
    if (arrow_angle == 270) {
      force = new Vec2(0,-1).times(100);
    }
    Vec2 hit_point = new Vec2(mouseX, mouseY);
    apply_force(force, hit_point, 1);
  }
  if (clicked_box2) {
    Vec2 force = new Vec2(1,0).times(100);
    if (arrow_angle == 0) {
      force = new Vec2(1,0).times(100);
    }
    if (arrow_angle == 90) {
      force = new Vec2(0,1).times(100);
    }
    if (arrow_angle == 180) {
      force = new Vec2(-1,0).times(100);
    }
    if (arrow_angle == 270) {
      force = new Vec2(0,-1).times(100);
    }
    Vec2 hit_point = new Vec2(mouseX, mouseY);
    apply_force(force, hit_point, 2);
  }
  
  ColideInfo info1 = collisionTest1(); //TODO: Use this result below
  
  //TODO the these values based on the results of a collision test
  Boolean hit_something1 = info1.hit; //Did I hit something?
  if (hit_something1){
    int hit_boxNum1 = info1.boxNum;
    Vec2 hit_point1 = info1.hitPoint;
    Vec2 hit_normal1 = info1.objectNormal;
    resolveCollision(hit_point1,hit_normal1,dt,hit_boxNum1);
  }
  
  ColideInfo info2 = collisionTest2(); //TODO: Use this result below
  
  //TODO the these values based on the results of a collision test
  Boolean hit_something2 = info2.hit; //Did I hit something?
  if (hit_something2){
    int hit_boxNum2 = info2.boxNum;
    Vec2 hit_point2 = info2.hitPoint;
    Vec2 hit_normal2 = info2.objectNormal;
    resolveCollision(hit_point2,hit_normal2,dt,hit_boxNum2);
  }
  
  ColideInfo objectInfo = objectCollisionTest(); //TODO: Use this result below
  
  //TODO the these values based on the results of a collision test
  Boolean object_hit_something = objectInfo.hit; //Did I hit something?
  if (object_hit_something){
    Vec2 hit_point = objectInfo.hitPoint;
    Vec2 hit_normal = objectInfo.objectNormal;
    resolveObjectCollision(hit_point,hit_normal,dt,objectInfo.boxNum);
  }
  
  background(200); //Grey background
  
  beginShape();
  texture(shelf);
  vertex(0, 0, 0, 0, 0);
  vertex(1200, 0, 0, shelf.width, 0);
  vertex(1200, 800, 0, shelf.width, shelf.height);
  vertex(0, 800, 0, 0, shelf.height);
  endShape();
  
  beginShape();
  texture(book);
  vertex(p3.x, p3.y, 0, 0, 0);
  vertex(p1.x, p1.y, 0, book.width, 0);
  vertex(p2.x, p2.y, 0, book.width, book.height);
  vertex(p4.x, p4.y, 0, 0, book.height);
  endShape();
  
  beginShape();
  texture(book);
  vertex(q3.x, q3.y, 0, 0, 0);
  vertex(q1.x, q1.y, 0, book.width, 0);
  vertex(q2.x, q2.y, 0, book.width, book.height);
  vertex(q4.x, q4.y, 0, 0, book.height);
  endShape();
  
  drawArrow(mouseX, mouseY, 100, arrow_angle);
}

void keyPressed(){
  if (key == 'r'){
    println("Resetting the simulation");
    momentum = new Vec2(0,0);          //Speed the box is translating (derivative of position)
    angular_momentum = 0;             //Speed the box is rotating (derivative of angle)

    center = new Vec2(400,400);        //Current position of center of mass
    angle = 0; /*radians*/          //Current rotation amount (orientation)
    
    momentum2 = new Vec2(0,0);          //Speed the box is translating (derivative of position)
    angular_momentum2 = 0;             //Speed the box is rotating (derivative of angle)

    center2 = new Vec2(800,400);        //Current position of center of mass
    angle2 = 0; /*radians*/          //Current rotation amount (orientation)
    
    return;
  }
  if (keyCode == UP) {
    arrow_angle = 270;
  }
  if (keyCode == LEFT) {
    arrow_angle = 180;
  }
  if (keyCode == RIGHT) {
    arrow_angle = 0;
  }
  if (keyCode == DOWN) {
    arrow_angle = 90;
  }
}

//Returns true iff the point 'point' is inside the box
boolean point_in_box(Vec2 point, Vec2 box_center, float box_w, float box_h, float box_angle){
  Vec2 relative_pos = point.minus(box_center);
  Vec2 box_right = new Vec2(cos(box_angle),sin(box_angle));
  Vec2 box_up = new Vec2(sin(box_angle),-cos(box_angle));
  float point_right = dot(relative_pos,box_right);
  float point_up = dot(relative_pos,box_up);
  if ((abs(point_right) < box_w/2) && (abs(point_up) < box_h/2))
    return true;
  return false;
}

void drawArrow(int cx, int cy, int len, float angle){
  pushMatrix();
  translate(cx, cy);
  rotate(radians(angle));
  line(-len,0,0, 0);
  line(0, 0,  - 8, -8);
  line(0, 0,  - 8, 8);
  popMatrix();
}

//---------------
//Vec 2 Library
//---------------

//Vector Library
//CSCI 5611 Vector 2 Library [Example]
// Stephen J. Guy <sjguy@umn.edu>

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
  
  public void normalize(){
    float magnitude = sqrt(x*x + y*y);
    x /= magnitude;
    y /= magnitude;
  }
  
  public Vec2 normalized(){
    float magnitude = sqrt(x*x + y*y);
    return new Vec2(x/magnitude, y/magnitude);
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

//2D cross product is a funny concept
// ...its the 3D cross product but with z = 0
// ... (only the resulting z compontent is not zero so we just store is as a scalar)
float cross(Vec2 a, Vec2 b){
  return a.x*b.y - a.y*b.x;
}

Vec2 projAB(Vec2 a, Vec2 b){
  return b.times(a.x*b.x + a.y*b.y);
}

Vec2 perpendicular(Vec2 a){
  return new Vec2(-a.y,a.x);
}
