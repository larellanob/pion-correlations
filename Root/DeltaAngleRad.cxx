Float_t DeltaAngle(Float_t x_rad, Float_t y_rad)
{
  Double_t x = x_rad*180./TMath::Pi();
  Double_t y = y_rad*180./TMath::Pi();

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
