//Solves the linear system Ax = b given A and b

float[] solveSystem(float[][] A, float[] b){
  int numIterations = 16;
  
  float[] x = b, newX = b;
  float transformedXcoord; //Helper variable for storing the relevant coordinate of (L + U)x
  
  for(int i = 0; i < numIterations; i++){
    
    // calculate new value for x = D^-1(b - (L + U)x)
    for(int j = 0; j < x.length; j++){
      
      //Calculate the relevant coordinate of L+U(x)
      transformedXcoord = 0;
      for(int k = 0; k < x.length; k++){
        if(k == j){
          break;
        }
        
        transformedXcoord += A[j][k] * x[k];
      }
      newX[j] = 1.0/A[j][j] * (b[j] - transformedXcoord);
    }
    if(i == numIterations - 1){
      float diff = 0.0;
      for(int k = 0; k < x.length; k++){
        float transformed = 0.0;
        for(int l = 0; l < x.length; l++){
          transformed += A[k][l]*newX[l];
        }
        diff += abs(b[k] - transformed);
      }
      println(diff);
    }
    x = newX;
  }
  
  return x;
}

//implements the jacobi method without having to construct the matrix A
float[] elementBasedSolver(float[] b, int sze){
  int iterationsPerCheckCycle = 1000;
  int maxCycles = 100;
  float[] x = new float[b.length], newX = new float[b.length];
  x = b.clone();
  float transformedXcoord; //Helper variable for storing the relevant coordinate of (L + U)x
  for(int checkCycle = 0; checkCycle < maxCycles; checkCycle++){
    for(int i = 0; i < iterationsPerCheckCycle; i++){
      // calculate new value for x = D^-1(b - (L + U)x)
      
      for(int j = 0; j < sze; j++){
        if(j == 0){
          transformedXcoord = x[j+1] + x[j + sze];
        }
        else if(j == sze-1){
          transformedXcoord = x[j-1] + x[j + sze];
        }
        else{
          transformedXcoord = x[j-1] + x[j+1] + x[j + sze];
        }
        
        newX[j] = 0.25 * (b[j] + transformedXcoord);
      }
      
      for(int j = sze; j < x.length-sze; j++){      
        if(j % sze == 0){
          transformedXcoord = x[j+1] + x[j - sze] + x[j + sze];
        }
        else if(j % sze == sze-1){
          transformedXcoord = x[j-1] + x[j - sze] + x[j + sze];
        }        
        else{
          transformedXcoord = x[j-1] + x[j+1] + x[j - sze] + x[j + sze];
        }
        
        newX[j] = 0.25 * (b[j] + transformedXcoord);
      }
      
      for(int j = sze * (sze-1); j < sze*sze; j++){
        if(j == sze * (sze - 1)){
          transformedXcoord = x[j+1] + x[j - sze];
        }
        else if(j == (sze*sze)-1){
          transformedXcoord = x[j-1] + x[j - sze];
        }
        else{
          transformedXcoord = x[j-1] + x[j+1] + x[j - sze];
        }
        
        newX[j] = 0.25 * (b[j] + transformedXcoord);
      }
      
      x = newX.clone();
    }
    
    if(checkAccuracy(b, x, sze) <= 0.5){
      println("Cycles: " + checkCycle);
      return x;
    }
  }
  println("Max Cycle count reached");
  println("Accuracy: " + checkAccuracy(b, x, sze));
  return x;
}

// Compute the euclidean distance between Ax and b
float checkAccuracy(float[] b, float[] x, int sze){
  float diff = 0.0;
  float approxBcoord; //Possible reuse of old values if all conditionals are false?
  for(int j = 0; j < sze; j++){
    if(j == 0){
      approxBcoord = -(x[j+1] + x[j + sze]) + 4 * x[j];
    }
    else if(j == sze-1){
      approxBcoord = -(x[j-1] + x[j + sze]) + 4 * x[j];
    }
    else{
      approxBcoord = -(x[j-1] + x[j+1] + x[j + sze]) + 4 * x[j];
    }
    
    diff += sq(b[j]-approxBcoord);
  }
  
  for(int j = sze; j < x.length-sze; j++){      
    if(j % sze == 0){
      approxBcoord = -(x[j+1] + x[j - sze] + x[j + sze]) + 4 * x[j];
    }
    else if(j % sze == sze-1){
      approxBcoord = -(x[j-1] + x[j - sze] + x[j + sze]) + 4 * x[j];
    }        
    else{
      approxBcoord = -(x[j-1] + x[j+1] + x[j - sze] + x[j + sze]) + 4 * x[j];
    }
    
    diff += sq(b[j]-approxBcoord);
  }
  
  for(int j = sze * (sze-1); j < sze*sze; j++){
    if(j == sze * (sze - 1)){
      approxBcoord = -(x[j+1] + x[j - sze]) + 4 * x[j];
    }
    else if(j == (sze*sze)-1){
      approxBcoord = -(x[j-1] + x[j - sze]) + 4 * x[j];
    }
    else{
      approxBcoord = -(x[j-1] + x[j+1] + x[j - sze]) + 4 * x[j];
    }
    diff += sq(b[j]-approxBcoord);
  }
  
  return sqrt(diff);
}
