simulation simul;
float tstep = 0.01;
int a = 0;
float[] zeroVec = {0, 0};

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
    du2[0][0] = (neighbors[1].du[0][0] - neighbors[3].du[0][0]) / 0.02;
    du2[0][1] = (neighbors[2].du[0][0] - neighbors[0].du[0][0]) / 0.02;
    
    du2[1][0] = (neighbors[1].du[1][1] - neighbors[3].du[1][1]) / 0.02;
    du2[1][1] = (neighbors[2].du[1][1] - neighbors[0].du[1][1]) / 0.02;
    
    //compute the diffusion term mu del^2 u
    float[] diffusion = {mu * (du2[0][0] + du2[0][1]), mu * (du2[1][0] + du2[1][1])};
      
    //compute du/dt according to incompressible navier-stokes equations ignoring pressure
    float[] dudt = new float[2];
    dudt[0] = diffusion[0] - convective[0];
    dudt[1] = diffusion[1] - convective[1];
    
    //if(abs(dudt[0]) > 100 || abs(dudt[1]) > 100){
    //  dudt[0] = 100 * abs(dudt[0])/dudt[0];
    //}
    
    //if(abs(dudt[1]) > 100){
    //  dudt[1] = 100 * abs(dudt[1])/dudt[1];
    //}
       
    if(magnitude(dudt) > 3){
      dudt[0] = 3 * dudt[0]/magnitude(dudt);
      dudt[1] = 3 * dudt[1]/magnitude(dudt);      
    }
        
    //float[] newVelocity = {u[0] + tstep*dudt[0], u[1] + tstep*dudt[1]};
    
    float[] newVelocity = {u[0] + dudt[0]*tstep, u[1] + dudt[1]*tstep}; //Not really new Velocity heheehheh
    if(magnitude(newVelocity) == 1.0/0.0){
      println("STOP!");
    }
    return new cell(x, y, newVelocity, mu);
  }
    
  void updatePartials(){
    du[0][0] = (neighbors[1].u[0] - neighbors[3].u[0]) / 0.02;
    du[0][1] = (neighbors[2].u[0] - neighbors[0].u[0]) / 0.02;
    
    du[1][0] = (neighbors[1].u[1] - neighbors[3].u[1]) / 0.02;
    du[1][1] = (neighbors[2].u[1] - neighbors[0].u[1]) / 0.02;
    
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
    float p = 0.0;
    
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        //vel[0] = 2*pow(2*pow(2,-1/(1-pow((i-150.0)/150,2)))*-sin(atan2(i-150, j-150)), 3);
        //vel[1] = 2*pow(2*pow(2,-1/(1-pow((j-150.0)/150,2)))*cos(atan2(i-150, j-150)), 3);
        
        float[] displacementFromCenter = {i-sze/2, j-sze/2};
        
        vel[1] = -sin(atan2(i-sze/2, j-sze/2)) * (1 - magnitude(displacementFromCenter)/(sze/2));
        vel[0] = cos(atan2(i-sze/2, j-sze/2)) * (1 - magnitude(displacementFromCenter)/(sze/2));
        cells[i][j] = new cell(i, j, vel, 0.001);
      }
    }
  }
  
  //void initNeighbors(){
  //  for(int i = 0; i < wth; i++){
  //    for(int j = 0; j < hgt; j++){
  //      cells[i][j].neighbors[0] = cells[i][(j+hgt-1) % hgt]; // up
  //      cells[i][j].neighbors[1] = cells[(i+wth+1) % wth][j]; // right
  //      cells[i][j].neighbors[2] = cells[i][(j+hgt+1) % hgt]; // down
  //      cells[i][j].neighbors[3] = cells[(i+wth-1) % wth][j]; // left
  //    }
  //  }
  //}
  
  void initNeighbors(){
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        if(j == 0){          
          cells[i][j].neighbors[0] = new cell(0, 0, zeroVec, 0);
        }
        else{
          cells[i][j].neighbors[0] = cells[i][j-1]; // up
        }

        if(j == sze-1){          
          cells[i][j].neighbors[2] = new cell(0, 0, zeroVec, 0);
        }
        else{
          cells[i][j].neighbors[2] = cells[i][j+1]; // down
        }
        
        if(i == 0){          
          cells[i][j].neighbors[3] = new cell(0, 0, zeroVec, 0);
        }
        else{
          cells[i][j].neighbors[3] = cells[i-1][j]; // left
        }

        if(i == sze-1){          
          cells[i][j].neighbors[1] = new cell(0, 0, zeroVec, 0);
        }
        else{
          cells[i][j].neighbors[1] = cells[i+1][j]; // left
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
    
    printArray(totalFlow);
    printArray(totalFlowMagnitude);
    
    pressureCorrection();
  }
  
  void pressureCorrection(){
    float[] velocityDivergenceField = new float[sze*sze]; //flattened list of divergences 
    
    //Find the rhs of the poisson equation
    
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        float duxdx = (cells[i][j].neighbors[1].u[0] - cells[i][j].neighbors[3].u[0]) / 0.02;
        float duydy = (cells[i][j].neighbors[2].u[0] - cells[i][j].neighbors[0].u[0]) / 0.02;
        velocityDivergenceField[i*sze + j] = duxdx + duydy;
      }
    }
    
    //Construct the Laplacian operator matrix
    
    float[][] laplacian = new float[sze*sze][sze*sze];
    for(int i = 0; i < sze*sze; i++){
      for(int j = 0; j < sze*sze; j++){
        if(i == j){
          laplacian[i][j] = 4;
        }
        
        else if(i == j + 1){
          laplacian[i][j] = -1;
        }
        
        else if(i == j - 1){
          laplacian[i][j] = -1;
        }
        
        else if(i == j + sze){
          laplacian[i][j] = -1;
        }
        
        else if(i == j - sze){
          laplacian[i][j] = -1;
        }
        
        else{
          laplacian[i][j] = 0;
        }
      }
    }
    
    float[] pField = solveSystem(laplacian, velocityDivergenceField);
    
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        float[] pgrad = new float[2];
        pgrad[0] = (pField[((i+1) * sze + j + sze*sze) % sze*sze] - pField[((i-1) * sze + j + sze*sze) % sze*sze])/ 0.02;
        pgrad[1] = (pField[(i * sze + j - 1 + sze*sze) % sze*sze] - pField[(i * sze + j + 1) % sze*sze])/ 0.02;
        
        cells[i][j].u[0] -= tstep*pgrad[0]; //eq. (6)
        cells[i][j].u[1] -= tstep*pgrad[1];
        //TODO: Implement eq. (6) to update velocities of cells
      }
    }
  }
  
  void display(){
    loadPixels();
    for(int i = 0; i < sze; i++){
      for(int j = 0; j < sze; j++){
        //println("coord: " +  i + "," + j);
        //println("vel: " + cells[i][j].v[0] + "," + cells[i][j].v[1]);
        pixels[j * sze + i] = HSVtoRGB(degrees(atan2(cells[i][j].u[0], cells[i][j].u[1])), magnitude(cells[i][j].u));
      }
    }
    updatePixels();
  }
}

void setup(){
  size(150, 150);
  simul = new simulation(150);
  simul.initNeighbors();
  print("Done!");
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
