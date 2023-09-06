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
    x = newX;
  }
  
  return x;
}
