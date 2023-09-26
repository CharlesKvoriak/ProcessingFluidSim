//DST algorithm as given in Chu, George 16.6

float[] DST(float[] input){
  float[][] x = new float[input.length*2+2][2];
  //create a vector x which is the input followed by itself flipped and negative. 0 inbetween them and at [0].
  for(int i = 0; i < input.length; i++){
    x[i+1][0] = input[i];
    x[2*input.length-i+1][0] = -input[i];
  }
  
  float[][] transformed = DFT(x);
  float[][] good = SFT(x);
  float[][] inverse = IDFT(transformed);
  
  float[] total = new float[2];
  for(int i = 0; i < good.length; i++){
    total[0] += good[i][0]*cos(-i*2.0*TWO_PI/(input.length));
    total[1] += good[i][1]*sin(-i*2.0*TWO_PI/(input.length));
  }
  
  println(good);
  println(inverse);
  float[] output = new float[input.length]; //<>//
  for(int i = 0; i < input.length; i++){
    output[i] = -0.5 * transformed[i + 1][1];
  }
  
  return output; //<>//
}

float[][] DFT(float[][] input){
  float[][] evens = new float[round(input.length/2.0)][2];
  float[][] odds = new float[round(input.length/2.0)][2];
  
  if(input.length == 1){
    return input;
  }
  
  for(int i = 0; i < input.length/2; i++){
    evens[i] = input[2*i];
    odds[i] = input[2*i + 1];
  }
  
  evens = DFT(evens);
  odds = DFT(odds);
  
  float[][] output = new float[input.length][2];
  
  for(int i = 0; i < input.length/2; i++){
    output[i] = add(evens[i], mult(euler(-i*TWO_PI/float(input.length)), odds[i]));
    output[i + input.length/2] = subtract(evens[i], mult(euler(-i*TWO_PI/float(input.length)), odds[i]));
  }

  return output;
}


float[] IDST(float[] input){
  float[][] x = new float[input.length*2+2][2];
  //create a vector x which is the input followed by itself flipped and negative. 0 inbetween them and at [0].
  for(int i = 0; i < input.length; i++){
    x[i+1][0] = input[i];
    x[2*input.length-i+1][0] = -input[i];
  }
  
  float[][] transformed = IDFT(x);

  float[] output = new float[input.length];
  for(int i = 0; i < input.length; i++){
    output[i] = 2.0/transformed[i + 1][1];
  }
  
  return output;
}

float[][] IDFT(float[][] input){
  float[][] evens = new float[round(input.length/2.0)][2];
  float[][] odds = new float[round(input.length/2.0)][2];
  
  if(input.length == 1){
    return input;
  }
  
  for(int i = 0; i < input.length/2; i++){
    evens[i] = input[2*i];
    odds[i] = input[2*i + 1];
  }
  
  evens = IDFT(evens);
  odds = IDFT(odds);
  
  float[][] output = new float[input.length][2];
  
  for(int i = 0; i < input.length/2; i++){
    float[] scaledown = {1.0/pow(float(input.length), i), 0};
    output[i] = add(evens[i], mult(scaledown, mult(euler(i*TWO_PI/float(input.length)), odds[i])));
    output[i + input.length/2] = subtract(evens[i], mult(scaledown, mult(euler(i*TWO_PI/float(input.length)), odds[i])));
  }

  return output;
}


float[] euler(float theta){
  float[] output = {cos(theta), sin(theta)};
  return output;
}

float[] mult(float[] a, float[] b){
  float[] output = {a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]};
  return output;
}

float[] add(float[] a, float[] b){
  float[] output = {a[0]+b[0], a[1]+b[1]};
  return output;
}

float[] subtract(float[] a, float[] b){
  float[] output = {a[0] - b[0], a[1] - b[1]};
  return output;
}
