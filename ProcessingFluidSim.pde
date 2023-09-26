//for vectors, first index is x and second is y
//for locations on the grid, first index is y and second is x

simulation simul;
float tstep = 0.01;
int a = 0;
float[] zeroVec = {0, 0};
float deltaX;

class cell{
  // clockwise [up, right, down, left]
  int x, y;
  float[] u = {0.0,0.0};
  float[][] du = {{0.0,0.0}, {0.0, 0.0}}; //partial derivatives of u
  float mu, p; //viscocity, pressure
  
  cell[] neighbors = new cell[4];
  cell(int xcoord, int ycoord, float[] velocity, float viscosity){
    x = xcoord;
    y = ycoord;
    arrayCopy(velocity, u);
    mu = viscosity;
    //printArray(velocity);
  }
  
  cell calculateProvisionalVelocity(){
    //compute the convective acceleration term (u \dot del)u
    float[] convective = {u[0] * du[0][0] + u[1] * du[1][0], u[0] * du[0][1] + u[1] * du[1][1]};
      
    //compute the second unmixed partial derivatives of u
    float[][] du2 = new float[2][2];
    du2[0][0] = (neighbors[1].du[0][0] - neighbors[3].du[0][0]) / (2* deltaX);
    du2[0][1] = (neighbors[2].du[0][0] - neighbors[0].du[0][0]) / (2* deltaX);
    
    du2[1][0] = (neighbors[1].du[1][1] - neighbors[3].du[1][1]) / (2* deltaX);
    du2[1][1] = (neighbors[2].du[1][1] - neighbors[0].du[1][1]) / (2* deltaX);
    
    //compute the diffusion term mu del^2 u
    float[] diffusion = {mu * (du2[0][0] + du2[0][1]), mu * (du2[1][0] + du2[1][1])};
      
    //compute du/dt according to incompressible navier-stokes equations ignoring pressure
    float[] dudt = new float[2];
    dudt[0] = diffusion[0] - convective[0];
    dudt[1] = diffusion[1] - convective[1];
       
    //if(magnitude(dudt) > 3){
    //  dudt[0] = 3 * dudt[0]/magnitude(dudt);
    //  dudt[1] = 3 * dudt[1]/magnitude(dudt);      
    //}
        
    //float[] newVelocity = {u[0] + tstep*dudt[0], u[1] + tstep*dudt[1]};
    
    float[] newVelocity = {u[0] + dudt[0]*tstep, u[1] + dudt[1]*tstep}; //Corresponds to u*
    return new cell(x, y, newVelocity, mu);
  }
    
  void updatePartials(){
    du[0][0] = (neighbors[1].u[0] - neighbors[3].u[0]) / (2 * deltaX);
    du[0][1] = (neighbors[2].u[0] - neighbors[0].u[0]) / (2 * deltaX);
    
    du[1][0] = (neighbors[1].u[1] - neighbors[3].u[1]) / (2 * deltaX);
    du[1][1] = (neighbors[2].u[1] - neighbors[0].u[1]) / (2 * deltaX);
    
    //dp[0] = (neighbors[1].p - neighbors[3].p) / 2;
    //dp[1] = (neighbors[2].p - neighbors[0].p) / 2;
  }
}

class simulation{
  cell[][] cells;
  int sze;
  
  simulation(int s){
    sze = s;
    float[] vel = {0, 0};
    cells = new cell[sze][sze];
    
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        //vel[0] = 2*pow(2*pow(2,-1/(1-pow((i-150.0)/150,2)))*-sin(atan2(i-150, j-150)), 3);
        //vel[1] = 2*pow(2*pow(2,-1/(1-pow((j-150.0)/150,2)))*cos(atan2(i-150, j-150)), 3);
        
        float[] displacementFromCenter = {i-sze/2, j-sze/2};
        
        vel[1] = -sin(atan2(j-sze/2, i-sze/2)) * exp(1 - 1/(1-sq(2*((magnitude(displacementFromCenter))/(sze/3.0))-1)));
        vel[0] = cos(atan2(j-sze/2, i-sze/2)) * exp(1 - 1/(1-sq(2*((magnitude(displacementFromCenter))/(sze/3.0))-1)));
        //vel[1] = 0.0;
        //vel[0] = -floor(2*j / sze);
        
        if(magnitude(displacementFromCenter) >= sze/3.0){
          vel[0] = 0;
          vel[1] = 0;
        }        
        
        if(vel[1] == 1.0/0.0){
         vel[1] = 0;
        }
        if(vel[0] == 1.0/0.0){
          vel[0] = 0;
        }
        
        cells[i][j] = new cell(i, j, vel, 0.0);
      }
    }
  }
  
  void initNeighbors(){
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        if(i == 0){          
          cells[i][j].neighbors[0] = cells[i][j];
        }
        else{
          cells[i][j].neighbors[0] = cells[i-1][j]; // up
        }

        if(i == sze-1){          
          cells[i][j].neighbors[2] = cells[i][j];
        }
        else{
          cells[i][j].neighbors[2] = cells[i+1][j]; // down
        }
        
        if(j == 0){          
          cells[i][j].neighbors[3] = cells[i][j];
        }
        else{
          cells[i][j].neighbors[3] = cells[i][j-1]; // left
        }

        if(j == sze-1){          
          cells[i][j].neighbors[1] = cells[i][j];
        }
        else{
          cells[i][j].neighbors[1] = cells[i][j+1]; // right
        }    
      }
    }
  }
  
  void update(){
    cell[][] newCells = new cell[sze][sze];
    float[] totalFlow = {0, 0};
    float totalFlowMagnitude = 0;
    
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        cells[i][j].updatePartials();
      }
    }
    
    for(int i = 0; i < sze; i++){
      arrayCopy(cells[i], newCells[i]);
      for(int j = 0; j < sze; j++){
        newCells[i][j] = cells[i][j].calculateProvisionalVelocity();
        totalFlow[0] += newCells[i][j].u[0];
        totalFlow[1] += newCells[i][j].u[1];
        totalFlowMagnitude += magnitude(newCells[i][j].u);
      }
    }
    
    for(int i = 0; i < sze; i++){
      arrayCopy(newCells[i] , cells[i]);
    }
    initNeighbors();
    
    //printArray(totalFlow);
    //printArray(totalFlowMagnitude);
    
    pressureCorrection();
  }
  
  void pressureCorrection(){
    float[] velocityDivergenceField = new float[sze*sze]; //flattened list of divergences 
    
    //Find the rhs of the poisson equation
    
    float total = 0;
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        //TODO: Fix boundaries by creating a special little function
        float duxdx = (cells[i][j].neighbors[1].u[0] - cells[i][j].neighbors[3].u[0]) / (2 * deltaX);
        float duydy = (cells[i][j].neighbors[2].u[1] - cells[i][j].neighbors[0].u[1]) / (2 * deltaX);
        velocityDivergenceField[i*sze + j] = -((4 * deltaX * deltaX)/tstep)*(duxdx + duydy);
        total += abs(duxdx + duydy);
      }
    }
    
    println("total divergence: " + total);
    
    //Construct the Laplacian operator matrix
    
    //float[][] laplacian = new float[sze*sze][sze*sze];
    //for(int i = 0; i < sze*sze; i++){
    //  for(int j = 0; j < sze*sze; j++){
    //    laplacian[i][j] = 0;
        
    //    if(i == j){
    //      laplacian[i][j] = -4;
    //    }
        
    //    if((j == i + 1) && (i % sze != sze-1)){
    //      laplacian[i][j] = 1;
    //    }
        
    //    if((j == i - 1) && (i % sze != 0)){
    //      laplacian[i][j] = 1;
    //    }
        
    //    if(i == j + sze){
    //      laplacian[i][j] = 1;
    //    }
        
    //    if(i == j - sze){
    //      laplacian[i][j] = 1;
    //    }
    //  }
    //}
    
    float[] pField = elementBasedSolver(velocityDivergenceField, sze).clone();
    
    
    println((cells[125][80+1].u[0] - cells[125][80-1].u[0] + cells[125+1][80].u[1] - cells[125-1][80].u[1])/(2*deltaX));
    
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        float[] pgrad = new float[2];
        // Boundaries are computed with forward/backward difference to avoid reaching out of the grid. All other values are computed with central difference.
        if(i == 0){
          pgrad[1] = (pField[((i+1) * sze) + j] - pField[(i*sze) + j])/deltaX;
        }
        
        else if(i == sze - 1){
          pgrad[1] = (pField[(i * sze) + j] - pField[((i-1)*sze) + j])/deltaX;
        }
        
        else{
          pgrad[1] = (pField[(i+1) * sze + j] - pField[(i-1) * sze + j])/ (2 * deltaX);
        }
        
        if(j == 0){
          pgrad[0] = (pField[i * sze + j + 1] - pField[i * sze + j])/deltaX;
        }
        
        else if(j == sze - 1){
          pgrad[0] = (pField[i * sze + j] - pField[i * sze + j - 1])/deltaX;
        }
        
        else{
          pgrad[0] = (pField[i * sze + j + 1] - pField[i * sze + j - 1])/ (2 * deltaX);
        }
        
        cells[i][j].u[0] -= (2*deltaX*tstep)*pgrad[0]; //eq. (6)
        cells[i][j].u[1] -= (2*deltaX*tstep)*pgrad[1];
      }
    }
    //int myX = 50;
    //int myY = 100;
    //println((cells[myY][myX+1].u[0] - cells[myY][myX-1].u[0])/(2*0.02) + (cells[myY+1][myX].u[1] - cells[myY-1][myX].u[1])/(2*0.02));
    
    total = 0;
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        float duxdx = (cells[i][j].neighbors[1].u[0] - cells[i][j].neighbors[3].u[0]) / (2 * deltaX);
        float duydy = (cells[i][j].neighbors[2].u[1] - cells[i][j].neighbors[0].u[1]) / (2 * deltaX);
        //if(i == 125 && j == 80){println(duxdx+duydy);}
        total += abs(duxdx + duydy);
      }
    }
    
    println("total divergence after correction: " + total);
    
    println();
    println((-4*pField[125*sze + 80] + pField[(125-1)*sze + 80] + pField[(125+1)*sze + 80] + pField[125*sze + 80 + 1] + pField[125*sze + 80 - 1])/(4 * deltaX * deltaX));
    //println(velocityDivergenceField[125*sze + 80]);
    println();
    println(velocityDivergenceField[125*sze + 80]/(-((4 * deltaX * deltaX))));
    println((cells[125][80+1].u[0] - cells[125][80-1].u[0] + cells[125+1][80].u[1] - cells[125-1][80].u[1])/(2*deltaX));
    println();
  }
  
  void display(){
    loadPixels();
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        //println("coord: " +  i + "," + j);
        //println("vel: " + cells[i][j].v[0] + "," + cells[i][j].v[1]);
        pixels[i * sze + j] = HSVtoRGB(degrees(atan2(cells[i][j].u[0], cells[i][j].u[1])), magnitude(cells[i][j].u));
      }
    }
    updatePixels();
  }
}

void setup(){
  size(150, 150);
  simul = new simulation(150);
  deltaX =  3.0 / float(150); //Constructed so that the simulation is always 3mx3m
  simul.initNeighbors();
  println("Done!");
  //float[] test = {0.0, 1.1, 2.2, 3.3, 4.4, 1.0, -4, 0.0};
  float[] test = {1, -2, 2};
  printArray(DST(test));
  printArray(IDST(DST(test)));
}

void draw(){
  simul.update();
  simul.display();
   
  a++;
  
  if(a % 100 == 0){
    print("hmmmmm");
  }
  //loadPixels();
  //for(int i = 0; i < width; i++){
  //  for(int j = 0; j < height; j++){
  //    pixels[j * width + i] = HSVtoRGB(degrees(atan2(i-150, j - 150)), sqrt(pow(i - 150, 2) + pow(j - 150, 2)) / 150);
  //  }
  //}
  //updatePixels();
}

float magnitude(float[] vector){
  return sqrt(pow(vector[0],2) + pow(vector[1],2));
}

//Not really HSV. My version (((:
color HSVtoRGB(float H, float S){
  //float H = radians(c[0]);
  //float S = c[1];
  H = radians(H);
  return color((cos(H) + 1)/ 2 * S * 255, (cos(H + TWO_PI/3) + 1)/2 * S * 255, (cos(H + 2 * TWO_PI/3) + 1)/2 * S * 255);
}
