Float_t DeltaAngle(Float_t x, Float_t y)
{

  Float_t result;

  result = abs(y-x);
  
  if ( result > 360. ) {
    cout << "Angle difference " << result << ", better check your code!" << endl;
  }
  if ( result > 180. ) {
    result = 360.-result;
  }

  return result;
}
