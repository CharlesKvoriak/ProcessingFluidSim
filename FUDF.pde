float[][] FFT(float[][] input){
  //println(input.length);
  if(input.length == 1){
    return input;
  }
  
  if((input.length & 1) == 1){ //if an odd number, to a standard Discrete Fourier Transform
    return SFT(input);
  }
  
  float[][] evens = new float[input.length/2][2];
  float[][] odds = new float[input.length/2][2];
  boolean parity = true; //true = even
  for(int i = 0; i<input.length; i++){
    if(parity){evens[i/2]=input[i];}
    else{odds[(i-1)/2]=input[i];}
    parity = !parity;
  }  
  
  float[][] TransformedEvens = FFT(evens);
  //println("evens: ");
  //for(int i = 0; i < TransformedEvens.length; i++){
  //  printArray(TransformedEvens[i]);
  //}
  float[][] TransformedOdds = FFT(odds);
  //println("odds: ");
  //for(int i = 0; i < TransformedOdds.length; i++){
  //  printArray(TransformedOdds[i]);
  //}  
  
  float[][] output = new float[input.length][2];
  parity = true;
  for(int i = 0; i<evens.length; i++){
    output[i] = add(TransformedEvens[i], mult(euler(TWO_PI/input.length*i), TransformedOdds[i]));
    output[i + (evens.length)] = add(TransformedEvens[i], mult(euler(TWO_PI/input.length*i+TWO_PI/2), TransformedOdds[i]));
  }
  
  return output; //<>//
}

float[][] IFT(float[][] input){
  float[][] output = new float [8][2];
  for(int i = 0; i<input.length; i++){
    float[] total = {0, 0};
    for(int j = 0; j<input.length; j++){
      total = add(total, mult(input[j], euler(TWO_PI*i*j/input.length)));
    }
    output[i] = total;
  }
  return output;
}

float[][] SFT(float[][] input){
  float[][] output = new float[input.length][2];
  for(int i = 0; i < input.length; i++){
    float[] total = {0, 0};
    for(int j = 0; j < input.length; j++){
      total = add(total, mult(input[j], euler(TWO_PI*i*j/input.length)));
    }
    output[i] = total;
  }
  
  return output;
}


//Inverse is almost identical to the FFT code

float[][] IFFT(float[][] input){
  //println(input.length);
  if(input.length == 1){
    return input;
  }
  
  if((input.length & 1) == 1){ //if an odd number, to a standard Inverse Discrete Fourier Transform
    return ISFT(input);
  }
  
  float[][] evens = new float[input.length/2][2];
  float[][] odds = new float[input.length/2][2];
  boolean parity = true; //true = even
  for(int i = 0; i<input.length; i++){
    if(parity){evens[i/2]=input[i];}
    else{odds[(i-1)/2]=input[i];}
    parity = !parity;
  }  
  
  float[][] TransformedEvens = IFFT(evens);
  //println("evens: ");
  //for(int i = 0; i < TransformedEvens.length; i++){
  //  printArray(TransformedEvens[i]);
  //}
  float[][] TransformedOdds = IFFT(odds);
  //println("odds: ");
  //for(int i = 0; i < TransformedOdds.length; i++){
  //  printArray(TransformedOdds[i]);
  //}  
  
  float[][] output = new float[input.length][2];
  parity = true;
  for(int i = 0; i<evens.length; i++){
    output[i] = add(TransformedEvens[i], mult(euler(-TWO_PI/input.length*i), TransformedOdds[i]));
    output[i + (evens.length)] = sub(TransformedEvens[i], mult(euler(-TWO_PI/input.length*i), TransformedOdds[i]));  
    //float[] scaledown = {1/float(input.length), 0};
    float[] scaledown = {1, 0};
    output[i] = mult(scaledown, output[i]);
    output[i + evens.length] = mult(scaledown, output[i + evens.length]);
  }
  
  return output;
}

float[][] ISFT(float[][] input){
  float[][] output = new float[input.length][2];
  for(int i = 0; i < input.length; i++){
    float[] total = {0, 0};
    for(int j = 0; j < input.length; j++){
      total = add(total, mult(input[j], euler(-TWO_PI*i*j/input.length)));
    }
    //float[] scaledown = {1/float(input.length), 0};
    float[] scaledown = {1, 0};
    output[i] = mult(total, scaledown);
    //output[i] =total;
  }
  
  return output;
}
