/* file:  ComputeTemperaturePrecipitationfrom018.cpp
 Scale present day monthly precipitation and temperature fields
 along the NGRIP oxygen isotope record.
 */

#include "./elements.h"
#include "../Numerics/numerics.h"

void ComputeD18OTemperaturePrecipitationFromPD(IssmDouble d018,IssmDouble dpermil,IssmDouble f,
					       IssmDouble* PrecipitationPresentday,IssmDouble* TemperaturePresentday,
					       IssmDouble* monthlytemperaturesout, IssmDouble* monthlyprecout){
  
  IssmDouble monthlytemperaturestmp[12],monthlyprectmp[12];
  IssmDouble deltaTemp;
  
  /* Constants */
  // dpermil = 2.4;/*degrees C per mil*/
  
  /*Create Delta Temp to be applied to monthly temps and used in precip scaling*/
  deltaTemp = dpermil * (d018+34.83);   
    
  for (int imonth = 0; imonth<12; imonth++){
    
    monthlytemperaturestmp[imonth] = TemperaturePresentday[imonth] + deltaTemp;
    monthlyprectmp[imonth] = PrecipitationPresentday[imonth]*exp((f/dpermil)*deltaTemp);
    
    /*Assign output pointer*/
    *(monthlytemperaturesout+imonth) = monthlytemperaturestmp[imonth];
    *(monthlyprecout+imonth) = monthlyprectmp[imonth];
  }
}
